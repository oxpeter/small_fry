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
import sys

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import Tkinter as tk   # for viewing the venn diagram
import Image, ImageDraw, ImageFont # for saving the venn diagram


from genomepy import config
import brain_machine as bm

########################################################################################

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "Performs set analysis on multiple lists and outputs venn diagrams and graphs showing how many are concordant/non-concordant etc")

    # input options
    parser.add_argument("experiments", metavar='filename', nargs='*', type=str,
                        help="a data file for comparing")
    parser.add_argument("orthologs", nargs=1, type=str,
                        help="orthomcl-formatted ortholog file")
    
    parser.add_argument("-m", "--mustcontain", type=str,
                        help="""comma-separated list containing 4-letter species codes
                        that must be present in each set of orthologs to include that set
                        """)
    parser.add_argument("-x", "--exclude", type=str,
                        help="""comma-separated list containing 4-letter species codes
                        that will be ignored in determining unique orthologous groups
                        """)
    parser.add_argument("-f", "--first", type=str,
                        help="specify the first DESeq2 P-value output file to use")
    parser.add_argument("-s", "--second", type=str,
                        help="specify the second DESeq2 P-value output file to use")
    parser.add_argument('-c', '--calibrate', type=str,
                        help="""provide an ortholog group name that is differentially
                        expressed in the same direction in all experiments, to ensure
                        that all directionality is consisent across experiments.""")
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="don't print messages")
    parser.add_argument("-o", "--output", type=str, default='lorf2.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-v", "--visible", action='store_true',default=False,
                        help="show charts and graphs")


    parser.add_argument("-L", "--list_genes", action='store_true',
                        help="""list common significant genes. (provide a file for name 
                        expansion""")
    parser.add_argument("-A", "--findall", action='store_true',
                        help="""find all significant concordant genes between given 
                        datasets""")
    parser.add_argument("-R", "--reversepolarity", action='store_true',
                        help="reverse the polarity of log fold change for second set")

    return parser



def fetch_orthologs(orthofile, mustcontain=None, exclude=None):
    handle = open(orthofile, 'rb')
    verbalise("M", "Converting from file", orthofile)
    ortho_dic = {}
    for line in handle:
        cols = line.split()
        if len(cols) > 0:
            counts = collections.Counter()
            for g in cols:
                spec = g.split('|')[0]
                if spec in exclude:
                    continue
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
                
    return ortho_dic

def translate_to_orthologs(degfile, orthodic):
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
                degdic[orthodic[cols[0]]] = padj, logfc
    handle.close()

    df = pd.DataFrame([ [k, v[0], v[1]] for k,v in degdic.items() ], 
                        columns=["gene","padj","logfc"])
    indexed_df = df.set_index('gene')
    return indexed_df

def replace_genenames(args):
    "Uses dictionary to replace and print only rows whose first column contains a match"
    ortho_dic = fetch_orthologs(args.orthologs)
    handle = open(args.second, 'rb')

    newfile = args.second[:-3] + 'ortho.out'
    while os.path.isfile(newfile): # make sure I'm not overwriting another file!
        newfile = newfile[:-4] + '1.out'
    orthofile = open( newfile, 'w')
    orthofile.write("baseMean log2FoldChange lfcSE stat pvalue padj\n")

    for line in handle:
        cols = line.split()
        if cols[0] in ortho_dic:
            cols[0] = ortho_dic[cols[0]]
            orthofile.write(" ".join(cols) + "\n")
    handle.close()
    orthofile.close()

    # update arg class so all downstream analyses use the newly created file:
    args.second = newfile

    return args

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

    verbalise("G", "Sum of circles", "1", len(pos1_u | neg1_u | concord_p | concord_n | discord_2p | discord_1p))
    verbalise("C", "Sum of circles", "2", len(pos2_u | neg2_u | concord_p | concord_n | discord_2p | discord_1p))


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
    verbalise("C", "Background concordance frequency = %.2f" % bkgd_freq)

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
    pp=PdfPages(outfile)
    plt.savefig(pp, format='pdf')
    pp.close()

    if visible:
        plt.show()
    else:
        plt.close()


def compare_two(args, logfile):
    df1 = pd.read_csv(args.first, sep=" ", header=0)
    df2 = pd.read_csv(args.second, sep=" ", header=0)
    df3 = df1.join(df2, lsuffix="1st").dropna()

    if len(df3.index)==0:
        verbalise("R", "No genes found in common between the two samples. If comparing different species, please make sure that you are using orthologous gene names")


    label1 = os.path.basename(args.first)[:5]
    label2 = os.path.basename(args.second)[:5]

    # get genes that are significant in both experiements and their direction:
    df1_sp = set(df3[(df3.padj1st<=0.05) & (df3.log2FoldChange1st>0)].index)
    df1_sn = set(df3[(df3.padj1st<=0.05) & (df3.log2FoldChange1st<0)].index)
    df1_nsp = set(df3[(df3.padj1st>0.05) & (df3.log2FoldChange1st>0)].index)
    df1_nsn = set(df3[(df3.padj1st>0.05) & (df3.log2FoldChange1st<0)].index)

    if not args.reversepolarity:
        df2_sp = set(df3[(df3.padj<=0.05) & (df3.log2FoldChange>0)].index)
        df2_sn = set(df3[(df3.padj<=0.05) & (df3.log2FoldChange<0)].index)
        df2_nsp = set(df3[(df3.padj>0.05) & (df3.log2FoldChange>0)].index)
        df2_nsn = set(df3[(df3.padj>0.05) & (df3.log2FoldChange<0)].index)
    else:
        df2_sp = set(df3[(df3.padj<=0.05) & (df3.log2FoldChange<0)].index)
        df2_sn = set(df3[(df3.padj<=0.05) & (df3.log2FoldChange>0)].index)
        df2_nsp = set(df3[(df3.padj>0.05) & (df3.log2FoldChange<0)].index)
        df2_nsn = set(df3[(df3.padj>0.05) & (df3.log2FoldChange>0)].index)


    print "DATAFRAMES:"
    verbalise("G", label1, "pos", len(df1_sp))
    verbalise("G", label1, "neg", len(df1_sn))
    verbalise("G", label1, "both", len(df1_sp | df1_sn))
    verbalise("C", label2, "pos", len(df2_sp))
    verbalise("C", label2, "neg", len(df2_sn))
    verbalise("C", label2, "both", len(df2_sp | df2_sn))

    print "CIRCLES"
    pos1_u,pos2_u,neg2_u,neg1_u,concord_p,discord_2p,concord_n,discord_1p = sigcounts(df1_sp, df2_sp, df1_sn, df2_sn)

    draw_circles(pos1_u,pos2_u,neg2_u,neg1_u,concord_p,discord_2p,concord_n,discord_1p,
                label1+" pos", label2+" pos", label1+" neg", label2+" neg",
                outfile=logfile[:-3] + "venn.png", visible=args.visible  )

    print "CONCORDANCE CHARTS"
    concordance_sets = concordancecounts(df1_sp, df2_sp, df1_sn, df2_sn,
                              df1_nsp, df2_nsp, df1_nsn, df2_nsn)


    draw_graph( *[len(s) for s in concordance_sets[:-1]],
                bkgd_freq=concordance_sets[-1], label1=label1, label2=label2,
                outfile=logfile[:-3] + "chart.pdf", visible=args.visible )

    return df1_sp & df2_sp, df1_sn & df2_sn

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

    if visible:
        root = tk.Tk()
        root.title('A Circle')

        # create canvases:
        canvas_1 = tk.Canvas(root, width=cw, height=ch, background='white')
        canvas_1.grid(row=0, column=1)
        # draw circles:
        canvas_1.create_oval(c1, outline='green', width=3)
        canvas_1.create_oval(c2, outline='red', width=3)
        canvas_1.create_oval(c3, outline='blue', width=3)
        canvas_1.create_oval(c4, outline='orange', width=3)

        # add text:
        text_c1 = canvas_1.create_text( 90,  110,  font=("Calibri", 14), anchor="nw", fill='darkgreen')
        text_c2 = canvas_1.create_text( 272, 110,  font=("Calibri", 14), anchor='nw', fill='red')
        text_c3 = canvas_1.create_text( 90,  292, font=("Calibri", 14), anchor='nw', fill='blue')
        text_c4 = canvas_1.create_text( 272, 292, font=("Calibri", 14), anchor='nw', fill='orange')

        text_o1 = canvas_1.create_text( 182, 130, font=("Calibri"), anchor='nw')
        text_o2 = canvas_1.create_text( 252, 200, font=("Calibri"), anchor='nw')
        text_o3 = canvas_1.create_text( 182, 272, font=("Calibri"), anchor='nw')
        text_o4 = canvas_1.create_text( 110, 200, font=("Calibri"), anchor='nw')

        canvas_1.itemconfig(text_c1, text=c1t)
        canvas_1.itemconfig(text_c2, text=c2t)
        canvas_1.itemconfig(text_c3, text=c3t)
        canvas_1.itemconfig(text_c4, text=c4t)
        canvas_1.itemconfig(text_o1, text=o1t)
        canvas_1.itemconfig(text_o2, text=o2t)
        canvas_1.itemconfig(text_o3, text=o3t)
        canvas_1.itemconfig(text_o4, text=o4t)

        # add labels:
        label_c1 = canvas_1.create_text( 70,  10,   font=("Calibri", 14), anchor="nw", fill='darkgreen')
        label_c2 = canvas_1.create_text( 242, 10,   font=("Calibri", 14), anchor="nw", fill='red')
        label_c3 = canvas_1.create_text( 242, 402,  font=("Calibri", 14), anchor="nw", fill='orange')
        label_c4 = canvas_1.create_text( 70,  402,  font=("Calibri", 14), anchor="nw", fill='blue')

        canvas_1.itemconfig(label_c1, text=l1t)
        canvas_1.itemconfig(label_c2, text=l2t)
        canvas_1.itemconfig(label_c3, text=l3t)
        canvas_1.itemconfig(label_c4, text=l4t)

        root.mainloop()

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

    image1.show()
    image1.save(outfile)
    image1.show()



########################################################################################

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    stop = False
    for f in args.experiments:
        if not os.path.isfile(f):
            verbalise("R", "%s could not be found" % (f))
            stop = True
    if stop:
        exit() 

    if args.exclude:
        exclusions = args.exclude.split(',')
    else:
        exclusions = None
        
    if args.mustcontain:
        necessary = args.mustcontain.split(',')
    else:
        necessary = None

    #args = replace_genenames(args)
    #verbalise("Y", "Ortholog file saved to %s" % os.path.basename(args.second))
    print "Ortholog file must contain:", args.mustcontain
    orthodic = fetch_orthologs(args.orthologs[0], 
                                mustcontain=necessary, 
                                exclude=exclusions)
    verbalise("G", "%d unique ortholog groups were found" % len(orthodic.items()))
    
    
    if args.first and args.second:
        # returns the genes that are significant and concordant in pos and neg directions:
        pcs, ncs = compare_two(args, logfile)
        common_to_all = pcs | ncs
    
    if args.experiments:

        # initialise containers:        
        concordant_array = pd.DataFrame(
                None,
                index=[os.path.basename(e)[:5] for e in args.experiments],
                columns=[os.path.basename(e)[:5] for e in args.experiments]
                )
        concordant_sig_genesets = []

        # do all pairwise comparisons:
        for (exp1, exp2) in itertools.product(args.experiments, args.experiments):
            # itertools.product produces an all-by-all comparison of the two lists
            # because this will include the same file compared to itself, we need
            # to remove those cases (which would mess up our comparisons):
            if exp1 == exp2:
                continue
            
            # get dataframes containing orthologs
            df1 = translate_to_orthologs(exp1, orthodic)
            df2 = translate_to_orthologs(exp2, orthodic)
            
            
            # the ratio of positive to negative log fold changes in a sample are 
            # probably consistent across species. So looking for a ratio gt or lt 1 in the
            # first sample, then convert all remaining samples to match, should give the
            # greatest concordance, and represent equivalent experimental conditions
            p1 = df1[(df1.logfc > 0) & (df1.padj <= 0.05)].count()["logfc"]
            n1 = df1[(df1.logfc < 0) & (df1.padj <= 0.05)].count()["logfc"]
            r1 = 1. * p1 / (n1 + 1)
            if args.calibrate and df1['logfc'].loc[args.calibrate] < 0:
                df1.logfc = df1.logfc * -1
            
            #if r1 < 1:
            #    df1.logfc = df1.logfc * -1
            
            p2 = df2[(df2.logfc > 0) & (df2.padj <= 0.05)].count()["logfc"]
            n2 = df2[(df2.logfc < 0) & (df2.padj <= 0.05)].count()["logfc"]
            r2 = 1. * p2 / (n2 + 1)
            if args.calibrate and df2['logfc'].loc[args.calibrate] < 0:
                df2.logfc = df2.logfc * -1
            #if r2 < 1:
            #    df1.logfc = df1.logfc * -1
            
                        
            df3 = df1.join(df2, how='left', lsuffix="1st").dropna()            
            label1 = os.path.basename(exp1)[:5]
            label2 = os.path.basename(exp2)[:5]
            


            # get genes that are significant in both experiements and their direction:
            df1_sp = set(df3[(df3.padj1st<=0.05) & (df3.logfc1st>0)].index)
            df2_sp = set(df3[(df3.padj   <=0.05) & (df3.logfc   >0)].index)
            df1_sn = set(df3[(df3.padj1st<=0.05) & (df3.logfc1st<0)].index)
            df2_sn = set(df3[(df3.padj   <=0.05) & (df3.logfc   <0)].index)

            df1_nsp = set(df3[(df3.padj1st>0.05) & (df3.logfc1st>0)].index)
            df2_nsp = set(df3[(df3.padj   >0.05) & (df3.logfc   >0)].index)
            df1_nsn = set(df3[(df3.padj1st>0.05) & (df3.logfc1st<0)].index)
            df2_nsn = set(df3[(df3.padj   >0.05) & (df3.logfc   <0)].index)

            verbalise(
    "B", 
    "\n\n%s (%d orthologs r %.1f) vs\n%s (%d orthologs r %.1f) : %d shared orthologs" % (
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
            verbalise("M", df3.loc['APR2016_00001325:'])
            concordance_sets = concordancecounts(df1_sp, df2_sp, df1_sn, df2_sn,
                              df1_nsp, df2_nsp, df1_nsn, df2_nsn)

            concordant_sig_genesets.append(concordance_sets[0])

            concordant_array[label1].loc[label2] = len(concordance_sets[0])
            verbalise("Y", "concordant DEGs = ", len(concordance_sets[0]))

        common_to_all = set.intersection(*concordant_sig_genesets)
        verbalise("R", "\nThere are %d genes common to all datasets" % len(common_to_all))

        print concordant_array
        
    if args.list_genes:
        """
        ncbi = bm.ncbi_dic(args.list_genes)
        genelist = open(logfile[:-3] + "concordantgenes.list", 'w')
        genelist.write("\n".join([bm.show_name(loc, ncbi) for loc in common_to_all]))
        genelist.close()
        print "\n".join([bm.show_name(loc, ncbi) for loc in common_to_all])
        """
        for o in common_to_all:
            verbalise("Y", o)