#!/usr/bin/env python

"""
This script compares the significant and non-significant genes between two experiments
and shows how many are concordant/non-concordant etc.
"""

import argparse
import collections
import itertools
import os
import re
import subprocess
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import Image, ImageDraw, ImageFont # for saving the venn diagram
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial import distance

from genomepy import config
from genomepy.genematch import Fisher_square, p_to_q
import brain_machine as bm

########################################################################################

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "Performs set analysis on multiple lists and outputs venn diagrams and graphs showing how many are concordant/non-concordant etc")

    # logging options
    parser.add_argument("-o", "--output", type=str, default='lorf2.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-D", "--display", action='store_true',default=False,
                        help="show charts and graphs")
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="don't print messages")

    # file input options
    parser.add_argument("experiments", metavar='filename', nargs='*', type=str,
                        help="a DESeq2 P-value output file")
    parser.add_argument("orthologs", nargs=1, type=str,
                        help="""orthomcl-formatted mcl output file. Each line of the file
                        starts with the ortholog group ID, followed by all members of the
                        group. The first 4 characters of each gene id identify the
                        species, followed by a pipe ("|") character. eg Cbir|LOC123456""")

    parser.add_argument("-m", "--mustcontain", type=str,
                        help="""comma-separated list containing 4-letter species codes
                        that must be present in each set of orthologs to include that set
                        """)
    parser.add_argument("-x", "--exclude", type=str,
                        help="""comma-separated list containing 4-letter species codes
                        that will be ignored in determining unique orthologous groups
                        """)
    parser.add_argument('-c', '--calibrate', type=str,
                        help="""provide an ortholog group name that is differentially
                        expressed in the same direction in all experiments, to ensure
                        that all directionality is consisent across experiments.""")

    # analysis options
    parser.add_argument("-L", "--list_genes", action='store_true',
                        help="""list all common significant genes.""")
    parser.add_argument('-n', '--name_genes', type=str,
                        help="""provide a file to convert gene loci to names.""")
    parser.add_argument("-A", "--findall", action='store_true',
                        help="""find all significant concordant genes between given
                        datasets""")
    parser.add_argument('-e', '--enrichment', type=str,
                        help="""provide a .obo file to perform GO term enrichment on the
                        common orthologs.""")


    return parser

def fetch_orthologs(orthofile, mustcontain=None, exclude=None):
    """
    Uses orthomcl's mcl output file to collect all appropriate ortholog groups. The
    dictionary created will only add groups when all members of the 'mustcontain' list
    are present, and when only a single representative of each present member exists,
    unless it is in the 'exclude' list, in which case any number of instances can be
    present.

    ortho_dic keys will be the gene name, and the values will be the ortholog name
    ortho_idx keys will be the ortholog name, and the values will be a list of all members

    NB: ortho_idx is not the strict subset that matches all necessary and exclusionary
    criteria! It is therefore not merely the inverse of ortho_dic. Instead, it contains
    all orthologs, allowing access to any of the ortholog names that were present in the
    original file.
    """
    handle = open(orthofile, 'rb')
    verbalise("M", "Converting from file", orthofile)
    ortho_dic = {}
    ortho_idx = {}
    for line in handle:
        cols = line.split()
        if len(cols) > 0:
            ortho_idx[cols[0]] = cols[1:]
            counts = collections.Counter()
            for g in cols:
                spec = g.split('|')[0]
                if spec in exclude:
                    counts[spec] = 1 # coerce to an acceptable number.
                else:
                    counts[spec] += 1
            if mustcontain:
                whitelist = mustcontain
            else:
                whitelist = counts.keys()

            # add only if all species present only once
            # and only if all necessary species are present
            if ( all(  c <= 1 for c in counts.values()) and
                 all( wl in counts.keys() for wl in whitelist )
               ):

                ortho_dic.update({ g:cols[0] for g in cols })
    handle.close()

    return ortho_dic, ortho_idx

def translate_to_orthologs(degfile, orthodic):
    """
    Takes the DESeq2 output file with genes, log2(fold change) and p-values, and
    creates a dictionary where the gene name is converted to the ortholog group name,
    to allow comparison between species.
    """
    degdic = {}
    handle = open(degfile, 'rb')
    for line in handle:
        cols = line.split()
        if len(cols) == 7 and cols[0] != "baseMean":
            if cols[0] in orthodic:
                try:
                    padj = float(cols[6])
                except ValueError:
                    if cols[5] == "NA":
                        padj = 1
                    elif float(cols[5]) < 0.05:
                        padj = 0.05
                    else:
                        padj = 1
                try:
                    logfc = float(cols[2])
                except ValueError:
                    logfc = 0

                try:
                    bmean = float(cols[1])
                except ValueError:
                    bmean = 0

                degdic[orthodic[cols[0]]] = padj, logfc, bmean
    handle.close()

    df = pd.DataFrame([ [k, v[0], v[1], v[2]] for k,v in degdic.items() ],
                        columns=["gene","padj","logfc","basemean"])
    indexed_df = df.set_index('gene')
    return indexed_df

def sigcounts(pos1, pos2, neg1, neg2):
    # get overlapping sets:
    concord_p = pos1 & pos2
    concord_n = neg1 & neg2
    discord_1p = pos1 & neg2
    discord_2p = pos2 & neg1

    # determine the genes unique to each set (non-overlapping):
    pos1_u = pos1 - concord_p - discord_1p
    pos2_u = pos2 - concord_p - discord_2p
    neg1_u = neg1 - concord_n - discord_2p
    neg2_u = neg2 - concord_n - discord_1p

    fordrawing = (pos1_u, pos2_u, neg2_u, neg1_u, concord_p, discord_2p, concord_n, discord_1p)
    return [str(len(s)) for s in fordrawing]

def concordancecounts(df1_sp, df2_sp, df1_sn, df2_sn, df1_nsp, df2_nsp, df1_nsn, df2_nsn):
    # signif in df 1:
    con_sig1   = (df1_sp & df2_sp)  | (df1_sn & df2_sn)
    con_nsig1  = (df1_sp & df2_nsp) | (df1_sn & df2_nsn)
    ncon_sig1  = (df1_sp & df2_sn)  | (df1_sn & df2_sp)
    ncon_nsig1 = (df1_sp & df2_nsn) | (df1_sn & df2_nsp)

    # signif in df 2:
    con_sig2   = (df1_sp & df2_sp)  | (df1_sn & df2_sn)
    con_nsig2  = (df1_nsp & df2_sp) | (df1_nsn & df2_sn)
    ncon_sig2  = (df1_sp & df2_sn)  | (df1_sn & df2_sp)
    ncon_nsig2 = (df1_nsp & df2_sn) | (df1_nsn & df2_sp)

    # background significance:
    bkgd_con = len((df1_nsp & df2_nsp) | (df1_nsn & df2_nsn) | (df1_sp & df2_sp) | (df1_sn & df2_sn))
    bkgd_dis = len((df1_nsp & df2_nsn) | (df1_nsn & df2_nsp) | (df1_sn & df2_sp) | (df1_sp & df2_sn))
    bkgd_freq = 1.0*bkgd_con/(bkgd_con+bkgd_dis)

    return (con_sig1, con_nsig1, ncon_sig1, ncon_nsig1,
            con_sig2, con_nsig2, ncon_sig2, ncon_nsig2,
            bkgd_freq)

def draw_graph( con_sig1, con_nsig1, ncon_sig1, ncon_nsig1,
                con_sig2, con_nsig2, ncon_sig2, ncon_nsig2,
                bkgd_freq=0.5, label1="group1", label2="group2",
                outfile="chart.pdf", visible=False ):

    N = 2
    ind = np.array([0.25,1.05])     # the x locations for the groups
    width = 0.5            # the width of the bars: can also be len(x) sequence

    totals = np.array([float(sum([con_sig1, con_nsig1, ncon_sig1, ncon_nsig1])),
                       float(sum([con_sig2, con_nsig2, ncon_sig2, ncon_nsig2]))])

    ra1 = np.array([con_nsig1, con_nsig2])
    ra2 = np.array([con_sig1, con_sig2])
    ra3 = np.array([ncon_nsig1, ncon_nsig2])
    ra4 = np.array([ncon_sig1, ncon_sig2])

    # define colors for chart:
    burgundy = (149./255,55./255,53./255)
    lightred = (217./255,150./255,148./255)
    greyblue = (85./255,142./255,213./255)
    lightblue = (142./255,180./255,227./255)

    linex = np.arange(0,1.9,0.3)
    liney = [ bkgd_freq for x in linex ]

    p1 = plt.bar(ind, ra1/totals, width, color=lightblue, alpha=1)
    p2 = plt.bar(ind, ra2/totals, width, color=greyblue, alpha=1, bottom=(ra1)/totals)
    p3 = plt.bar(ind, ra3/totals, width, color=lightred, alpha=1, bottom=(ra1+ra2)/totals)
    p4 = plt.bar(ind, ra4/totals, width, color=burgundy, alpha=1, bottom=(ra1+ra2+ra3)/totals)
    plt.plot(linex, liney, color='black', linewidth=2, linestyle='--')

    plt.ylabel('frequency')
    plt.title('concordance of gene expression data between experiments')
    plt.xticks(ind+width/2., ["%s\n%d DEGs" % (label1, totals[0]), "%s\n%d DEGs" % (label2,totals[1])])
    plt.yticks(np.arange(0,1,.20))

    def autolabel(rects, heights, adjs, values):
        # attach some text labels
        for i, rect in enumerate(rects):
            plt.text(rect.get_x()+rect.get_width()/2,
                    heights[i]-adjs[i]/2-0.03,
                    '%d'%(int(values[i])),
                    ha='center', va='bottom', size=16, weight='bold')

    autolabel(p1, (ra1)/totals,             ra1/totals, ra1)
    autolabel(p2, (ra1+ra2)/totals,         ra2/totals, ra2)
    autolabel(p3, (ra1+ra2+ra3)/totals,     ra3/totals, ra3)
    autolabel(p4, (ra1+ra2+ra3+ra4)/totals, ra4/totals, ra4)

    plt.savefig(outfile, format='pdf')


    if visible:
        plt.show()
    else:
        plt.close()

def draw_circles(c1t,c2t,c3t,c4t,o1t,o2t,o3t,o4t,l1t,l2t,l3t,l4t, outfile="venn.jpg", visible=False):
    # define colors:
    burgundy =  (149,55,53)
    lightred =  (217,150,148)
    greyblue =  (85,142,213)
    lightblue = (142,180,227)
    green =     (5,128,0)
    orange =    (228,108,9)
    black =     (0,0,0)
    white =     (255,255,255)

    # set canvas size:
    cw = 425
    ch = 480

    # specify top left and bottom right coords of circles:
    c1 = (20,40,220,240)
    c2 = (162,40,362,240)
    c3 = (20,182,220,382)
    c4 = (162,182,362,382)

    # create canvases:
    image1 = Image.new("RGB", (cw, ch), white)
    draw = ImageDraw.Draw(image1)

    # draw circles:
    draw.ellipse(c1, fill=None, outline=green)
    draw.ellipse(c2, fill=None, outline=burgundy)
    draw.ellipse(c3, fill=None, outline=greyblue)
    draw.ellipse(c4, fill=None, outline=orange)


    # add text:
    draw.text((90,  110),c1t,green, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 18))
    draw.text((272, 110),c2t,burgundy, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 18))
    draw.text((90,  292),c3t,greyblue, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 18))
    draw.text((272, 292),c4t,orange, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 18))

    draw.text((182, 130),o1t,black, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 18))
    draw.text((252, 200),o2t,black, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 18))
    draw.text((182, 272),o3t,black, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 18))
    draw.text((110, 200),o4t,black, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 18))

    # add labels:
    draw.text((70,  10),l1t,green, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 16))
    draw.text((242, 10),l2t,burgundy, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 16))
    draw.text((242, 402),l3t,orange, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 16))
    draw.text((70,  402),l4t,greyblue, font=ImageFont.truetype('/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif-Bold.ttf', 16))

    image1.save(outfile)
    if visible:
        image1.show()

    return outfile

def distance_tree(array, outfile, metric='euclidean', method='complete',
                    basis="shared degs", display=True):

    verbalise("C", "%s" % (basis))
    verbalise(str(array))
    print "\n"

    inverse = 1 / array
    zeroed  = inverse.fillna(0)
    Y  = distance.pdist(zeroed, metric)
    Z  = hierarchy.linkage(Y, method=method, metric=metric)

    plt.figure(figsize=(10,7))
    hierarchy.dendrogram(Z,
                        leaf_rotation=90.,
                        leaf_font_size=8.,
                        labels=array.columns)
    plt.title(
        "Distances based on %s \n(%s distance tree, %s clustering)" % (
                        basis, metric, method)
              )
    plt.ylabel("Distance")
    plt.savefig(outfile, format='pdf')
    if display:
        plt.show()
    else:
        plt.close()

def pairwise_container(names, truncate=True):
    if truncate:
        df = pd.DataFrame(None,
                    index=[os.path.basename(e)[:5] for e in names],
                    columns=[os.path.basename(e)[:5] for e in names]
                    )
    else:
        df = pd.DataFrame(None,
                    index=[os.path.basename(e) for e in names],
                    columns=[os.path.basename(e) for e in names]
                    )

    return df

########################################################################################

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    # check all input files can be found
    stop = False
    for f in args.experiments:
        if not os.path.isfile(f):
            verbalise("R", "%s could not be found" % (f))
            stop = True
    if stop:
        sys.exit(1)

    # set up lists of species that must be included or can be ignored in ortholog groups
    if args.exclude:
        exclusions = args.exclude.split(',')
    else:
        exclusions = None

    if args.mustcontain:
        necessary = args.mustcontain.split(',')
    else:
        necessary = None

    verbalise("B", "\nFinding all orthologs containing:", " ".join(args.mustcontain))
    verbalise("B", "and excluding:", " ".join(args.exclude))
    orthodic, ortho_idx = fetch_orthologs(args.orthologs[0],
                                            mustcontain=necessary,
                                            exclude=exclusions)
    verbalise("G", "%d genes were indexed in %d ortholog groups" % (len(orthodic),
                                                                    len(ortho_idx) ))

    if len(args.experiments) > 2:

        ########################################################################
        ######## combine all experiments into single dataframe #################
        dfall = translate_to_orthologs(args.experiments[0], orthodic)
        for num, exp in enumerate(args.experiments[1:]):
            dfnew = translate_to_orthologs(exp, orthodic)
            dfall = dfall.join(dfnew, how='left', rsuffix="_%d" % (num+1)).dropna()

        verbalise("C", "%d orthologous groups meet selection criteria" % len(dfall))
        verbalise("M", "%d orthologs with >= 2 fold difference in at least 1 of %d species" %
                 (len(dfall[(dfall.logfc   >= 1.0) |(dfall.logfc   <= -1.0) | \
                            (dfall.logfc_1 >= 1.0) |(dfall.logfc_1 <= -1.0) | \
                            (dfall.logfc_2 >= 1.0) |(dfall.logfc_2 <= -1.0) | \
                            (dfall.logfc_3 >= 1.0) |(dfall.logfc_3 <= -1.0) | \
                            (dfall.logfc_4 >= 1.0) |(dfall.logfc_4 <= -1.0) | \
                            (dfall.logfc_5 >= 1.0) |(dfall.logfc_5 <= -1.0) ]),
                  len(args.experiments)) )
        verbalise("M", "%d orthologs with negative >= 2 fold difference in %d species" %
                 (len(dfall[(dfall.logfc   < 0) & \
                            (dfall.logfc_1 < 0) & \
                            (dfall.logfc_2 < 0) & \
                            (dfall.logfc_3 < 0) & \
                            (dfall.logfc_4 < 0) & \
                            (dfall.logfc_5 < 0) ]),
                  len(args.experiments)) )
        verbalise("M", "%d orthologs with positive >= 2 fold difference in %d species" %
                 (len(dfall[(dfall.logfc   > 0)  & \
                            (dfall.logfc_1 > 0)  & \
                            (dfall.logfc_2 > 0)  & \
                            (dfall.logfc_3 > 0)  & \
                            (dfall.logfc_4 > 0)  & \
                            (dfall.logfc_5 > 0)  ]),
                  len(args.experiments)) )

        names =  [ os.path.basename(e)[:5] for e in args.experiments ]
        cols = [ 'logfc' ] + [ 'logfc_%d' % e for e in range(1,len(args.experiments)) ]
        converter = dict(zip(cols, names))

        special_array = pairwise_container(names, truncate=False)
        for (exp1, exp2) in itertools.product(cols, repeat=2):
            if exp1 == exp2:
                continue
            corr = dfall[exp1].corr(dfall[exp2], method='pearson')
            special_array[converter[exp1]].loc[converter[exp2]] = corr

        # calculate distance based on number of shared DEGs:
        distance_tree(special_array,
                        "%sspecial_tree.pdf" % logfile[:-3],
                        metric='euclidean',
                        method='single',
                        basis="all genes with fold change >= 2",
                        display=args.display)

        #############################################
        ######## do all pairwise comparisons ########

        ######## initialise containers ##############
        # for storing pairwise shared DEGs:
        concordant_array = pairwise_container(args.experiments)
        # for storing pairwise correlation of log2(fold change):
        correl_array = pairwise_container(args.experiments)
        # for storing (inverse) jaccard's index of DEGs:
        jaccard_array = pairwise_container(args.experiments)
        # for storing correlation of highly differential genes:
        high_correl_array = pairwise_container(args.experiments)
        # for storing binary correlation of all genes:
        bit_correl_array = pairwise_container(args.experiments)
        # for storing number of common highly differential genes:
        diff_count_array = pairwise_container(args.experiments)
        # for storing number of concordant genes
        all_conc_array = pairwise_container(args.experiments)
        # for storing number of concordant genes significant in at least one species
        goodenough_array = pairwise_container(args.experiments)

        concordant_sig_genesets = []

        # for storing graphs and charts:
        pdfhandle = PdfPages(logfile[:-3] + "barcharts.pdf")
        all_pngs = []

        for (exp1, exp2) in itertools.product(args.experiments, args.experiments):
            """
            itertools.product produces an all-by-all comparison of the two lists
            because this will include the same file compared to itself, we need
            to remove those cases (which would mess up our comparisons):
            """
            if exp1 == exp2:
                continue

            # get dataframes containing orthologs
            df1 = translate_to_orthologs(exp1, orthodic)
            df2 = translate_to_orthologs(exp2, orthodic)

            # if ortholog group is provided for polarity calibration, make sure all
            # changes for this ortholog are positive.
            if args.calibrate and df1['logfc'].loc[args.calibrate] < 0:
                df1.logfc = df1.logfc * -1
            if args.calibrate and df2['logfc'].loc[args.calibrate] < 0:
                df2.logfc = df2.logfc * -1

            # calculate the ratio of positive to negative changes
            p1 = df1[(df1.logfc > 0) & (df1.padj <= 0.05)].count()["logfc"]
            n1 = df1[(df1.logfc < 0) & (df1.padj <= 0.05)].count()["logfc"]
            r1 = 1. * p1 / (n1 + 1)

            p2 = df2[(df2.logfc > 0) & (df2.padj <= 0.05)].count()["logfc"]
            n2 = df2[(df2.logfc < 0) & (df2.padj <= 0.05)].count()["logfc"]
            r2 = 1. * p2 / (n2 + 1)

            # create merged dataset with only genes common to both species:
            df3 = df1.join(df2, how='left', lsuffix="1st").dropna()
            label1 = os.path.basename(exp1)[:5]
            label2 = os.path.basename(exp2)[:5]

            ######### analyse pairwise gene sets #######
            # calculate various correlation metrics:
            correlation = df3['logfc'].corr(df3['logfc1st'], method='pearson')

            high_expr = 1  #0.84799690655495  =log2(1.8)
            high_corr   = df3.loc[                                                 \
                            ((df3.logfc >= high_expr) | (df3.logfc <= -high_expr)) &     \
                            ((df3.logfc1st >= high_expr) | (df3.logfc1st <= -high_expr)) \
                             ]['logfc'].corr(df3.loc[                              \
                                                    ((df3.logfc >= high_expr) |
                                                    (df3.logfc <= -high_expr) ) &  \
                                                    ((df3.logfc1st >= high_expr) |
                                                    (df3.logfc1st <= -high_expr)) \
                                                     ]['logfc1st'], method='pearson')

            hsize = df3.loc[((df3.logfc >= high_expr)    | (df3.logfc <= -high_expr)) &  \
                            ((df3.logfc1st >= high_expr) | (df3.logfc1st <= -high_expr)) \
                            ].count()['logfc']


            df3['bitfc'] = df3.apply(lambda x: x['logfc']/abs(x['logfc']), axis=1)
            df3['bitfc1st'] = df3.apply(lambda x: x['logfc1st']/abs(x['logfc1st']), axis=1)

            bit_corr    = df3['bitfc'].corr(df3['bitfc1st'], method='pearson')

            # store correlation values in containers
            correl_array[label1].loc[label2]        = correlation
            high_correl_array[label1].loc[label2]   = high_corr
            bit_correl_array[label1].loc[label2]    = bit_corr
            diff_count_array[label1].loc[label2]    = hsize

            # get genes that are significant in both experiements and their direction:
            df1_sp = set(df3[(df3.padj1st<=0.05) & (df3.logfc1st>0)].index)
            df2_sp = set(df3[(df3.padj   <=0.05) & (df3.logfc   >0)].index)
            df1_sn = set(df3[(df3.padj1st<=0.05) & (df3.logfc1st<0)].index)
            df2_sn = set(df3[(df3.padj   <=0.05) & (df3.logfc   <0)].index)

            df1_nsp = set(df3[(df3.padj1st>0.05) & (df3.logfc1st>0)].index)
            df2_nsp = set(df3[(df3.padj   >0.05) & (df3.logfc   >0)].index)
            df1_nsn = set(df3[(df3.padj1st>0.05) & (df3.logfc1st<0)].index)
            df2_nsn = set(df3[(df3.padj   >0.05) & (df3.logfc   <0)].index)

            # get number of genes concordantly expressed (significant or not)
            num_concord = len(df3[((df3.logfc1st > 0) & (df3.logfc > 0)) |
                                  ((df3.logfc1st < 0) & (df3.logfc < 0))])
            all_conc_array[label1].loc[label2] = num_concord

            # get number of genes concordantly expressed (significant in at least one)
            num_close = len(df3[((df3.logfc1st > 0) & (df3.logfc > 0) |
                                (df3.logfc1st > 0) & (df3.logfc > 0)) &
                                ((df3.padj <= 0.05)  |  (df3.padj1st <= 0.05))])
            goodenough_array[label1].loc[label2] = num_close

            # find numbers of concordant and discordant genes:
            concordance_sets = concordancecounts(df1_sp, df2_sp, df1_sn, df2_sn,
                              df1_nsp, df2_nsp, df1_nsn, df2_nsn)
            concordant_sig_genesets.append(concordance_sets[0])
            concordant_array[label1].loc[label2] = len(concordance_sets[0])

            # find inverse of jaccard's index for significant genes:
            ijidx = 1. *(len(df1_sp | df1_sn) + len(df2_sp | df2_sn) - len(concordance_sets[0])) / \
                    (len(df1_sp | df1_sn) + len(df2_sp | df2_sn) - 2 * len(concordance_sets[0]))
            jaccard_array[label1].loc[label2] = ijidx


            ########## report output ###################
            verbalise("B",
    "\n\n%s (%d orthologs, er=%.1f) vs\n%s (%d orthologs, er=%.1f) : %d shared orthologs" % (
                     label1.upper(), len(df1), r1,
                     label2.upper(), len(df2), r2,
                     len(df3))
                      )
            verbalise("G", "%s: %d significant and pos" % (label1, len(df1_sp)))
            verbalise("G", "%s: %d significant and neg" % (label1, len(df1_sn)))
            verbalise("G", "%s: %d all significant    " % (label1, len(df1_sp | df1_sn)))
            verbalise("C", "%s: %d significant and pos" % (label2, len(df2_sp)))
            verbalise("C", "%s: %d significant and neg" % (label2, len(df2_sn)))
            verbalise("C", "%s: %d all significant    " % (label2, len(df2_sp | df2_sn)))
            verbalise("Y", "concordant DEGs = ", len(concordance_sets[0]))
            verbalise("M", "%d genes w log2fc > %d" % (hsize, high_expr))

            # create graphs of overlapping DEGs:
            scounts = sigcounts(df1_sp, df2_sp, df1_sn, df2_sn)
            pos1_u,pos2_u,neg2_u,neg1_u,concord_p,discord_2p,concord_n,discord_1p = scounts

            outpng = draw_circles(pos1_u, pos2_u,
                        neg2_u, neg1_u,
                        concord_p, discord_2p,
                        concord_n, discord_1p,
                        label1+" pos", label2+" pos", label1+" neg", label2+" neg",
                        outfile= "%s%s_%s.venn.png" % (logfile[:-3],label1,label2 ),
                        visible=args.display  )

            all_pngs.append(outpng) # will be concatenated into a single pdf later.
            concordance_sets = concordancecounts(df1_sp, df2_sp, df1_sn, df2_sn,
                                                 df1_nsp, df2_nsp, df1_nsn, df2_nsn)
            draw_graph( *[len(s) for s in concordance_sets[:-1]],
                        bkgd_freq=concordance_sets[-1],
                        label1=label1, label2=label2,
                        outfile=pdfhandle,
                        visible=args.display )

        # collate images, close output files and cleanup:
        pdfhandle.close()

        # convert is from the ImageMagick ubuntu install
        status = subprocess.call(
                ["convert"] + all_pngs + ["%svenn_diagrams.pdf" % logfile[:-3]]
                                )
        for f in all_pngs:
            os.remove(f)

        ########### SUMMARY OF ALL DATASETS ###########
        common_to_all = set.intersection(*concordant_sig_genesets)

        # calculate distance based on number of shared DEGs:
        distance_tree(concordant_array,
                        "%sdeg_tree.pdf" % logfile[:-3],
                        metric='euclidean',
                        method='single',
                        basis="number of concordant DEGs",
                        display=args.display)

        # calculate distance based on number of all concordant genes:
        distance_tree(all_conc_array,
                        "%sall_concordant_tree.pdf" % logfile[:-3],
                        metric='euclidean',
                        method='single',
                        basis="number of concordant genes",
                        display=args.display)

        # calculate distance based on number of all concordant genes with >= 1 sig:
        distance_tree(goodenough_array,
                        "%sall_concordant_gte1_sig_tree.pdf" % logfile[:-3],
                        metric='euclidean',
                        method='single',
                        basis="number of concordant genes, at least 1 significant",
                        display=args.display)



        # calculate distance based on jaccard's index of DEGs:
        distance_tree(jaccard_array,
                        "%sjaccards_tree.pdf" % logfile[:-3],
                        metric='euclidean',
                        method='single',
                        basis="Jaccard's index of DEGs",
                        display=args.display)

        # calculate distance based on log2(fold change)
        distance_tree(correl_array,
                        "%slog2fc_correl_tree.pdf" % logfile[:-3],
                        metric='euclidean',
                        method='single',
                        basis="correlation of log2(fold change)",
                        display=args.display)

        # calculate distance based on log2(fold change) of highly expressed genes
        distance_tree(high_correl_array,
                        "%slog2fc_high_correl_tree.pdf" % logfile[:-3],
                        metric='euclidean',
                        method='single',
                        basis="correlation of log2(fold change) of highly differential genes",
                        display=args.display)

        # calculate distance based on bitwise log2(fold change)
        distance_tree(bit_correl_array,
                        "%slog2fc_bit_correl_tree.pdf" % logfile[:-3],
                        metric='euclidean',
                        method='single',
                        basis="correlation of bitwise log2(fold change)",
                        display=args.display)

        # calculate distance based on bitwise log2(fold change)
        distance_tree(diff_count_array,
                        "%slog2fc_highdiff_count_tree.pdf" % logfile[:-3],
                        metric='euclidean',
                        method='single',
                        basis="number of common highly differential genes",
                        display=args.display)




        verbalise("R", "\nThere are %d genes common to all datasets" % len(common_to_all))

    else:
        verbalise("R", "Insufficient output files were provided")
        sys.exit(1)


    # list all orthologs that were concordantly expressed in all pairwise comparisons:
    if args.list_genes:
        if args.name_genes:     # create dictionary for finding gene names
            name_chart = {}
            handle = open(args.name_genes, 'rb')
            for line in handle:
                cols = line.split()
                name_chart[cols[0]] = " ".join(cols[1:])
            handle.close()


        for o in common_to_all:
            name = "---"
            if args.name_genes:
                if o in ortho_idx:
                    for gene in ortho_idx[o]:
                        if gene[5:] in name_chart:
                            name = name_chart[gene[5:]]

            if name == "---":
                name = "Includes: %s" % " ".join(ortho_idx[o][:3])
            verbalise("Y", "%20s %s" % (o, name))


    # perform Fisher's Exact Test for enrichment of Pfam domains in the orthologs that
    # were concordantly expressed in all pairwise comparisons
    if args.enrichment:
        background_set_size = len(df3)
        background_gene_set = { g[5:]:True for o in df3.index for g in ortho_idx[o]}
        background_gene_set.update({ g:True for o in df3.index for g in ortho_idx[o]})


        # create pfam indices:
        pfam_idx = {}   # store all genes that contain a given pfam acc.
        pfam_ref = {}   # store definitions of pfam accession numbers
        pfam_acc = {}   # store all pfams for a given gene

        handle = open(args.enrichment, 'rb')
        for line in handle:
            cols = line.split()

            # skip all null lines:
            if len(cols) == 0:
                continue
            elif line[0] == '#':
                continue

            # add to all the indices
            if cols[1] not in pfam_ref:
                pfam_ref[cols[1]] = " ".join(cols[22:])
                pfam_idx[cols[1]] = []

            if cols[3] not in pfam_acc:
                pfam_acc[cols[3]] = []

            pfam_idx[cols[1]].append(cols[3])
            pfam_acc[cols[3]].append(cols[1])

        handle.close()

        ### create fisher squares for all concordant common orthologs:

        # create key sets and counters:
        common_pfams = set()
        fs = [] # contains all the fishers tests

        # collect all pfams in the CCOs:
        deg_pfam_counter = collections.Counter()
        for o in common_to_all:
            for gene in ortho_idx[o]:
                if gene[5:] in pfam_acc:
                    common_pfams.update(pfam_acc[gene[5:]])
                    deg_pfam_counter.update(set(pfam_acc[gene[5:]]))
                    break
                elif gene in pfam_acc:
                    common_pfams.update(pfam_acc[gene])
                    deg_pfam_counter.update(set(pfam_acc[gene]))
                    break

        # count occurrences of pfam in each category
        for pfam in common_pfams:
            degs_w_pfam  = deg_pfam_counter[pfam]
            degs_wo_pfam = len(common_to_all) - degs_w_pfam
            nogs_w_pfam  = sum( 1 for g in pfam_idx[pfam] if g in background_gene_set) - degs_w_pfam
            nogs_wo_pfam = background_set_size - len(common_to_all) - nogs_w_pfam

            fs.append(Fisher_square(TwG=degs_w_pfam,
                                    TwoG=degs_wo_pfam,
                                    NwG=nogs_w_pfam,
                                    NwoG=nogs_wo_pfam,
                                    id=pfam,
                                    id_info="%s (%s)" % (pfam, pfam_ref[pfam]) )
                    )

        for s in fs:
            s.fisher_test()
        pvals = [ s.p for s in fs ]
        fdr = p_to_q(pvals, alpha=0.05, method='fdr_bh')
        verbalise("M", "%d pfam domains found in concordant common orthologs" % (len(fs)))

        for s in fs:
            if s.p in fdr:
                if fdr[s.p] <= 0.05 and s.TwG > 1:
                    verbalise("C", s)
                    verbalise("M", "FDR = %.3f" % fdr[s.p])
            elif s.p <= 0.05 and s.TwG > 1:
                verbalise("C", s)
                verbalise("M", "FDR not determined for P-value %r" % s.p)
