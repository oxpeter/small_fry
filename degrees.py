#!/usr/bin/env python

"""
This script compares the significant and non-significant genes between two experiments
and shows how many are concordant/non-concordant etc.
"""
import os
import sys
import re
import argparse
import itertools

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import Tkinter as tk

from genomepy import config
import brain_machine as bm

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "Performs set analysis on multiple lists and outputs venn diagrams.")

    # input options
    parser.add_argument("experiments", metavar='filename', nargs='*', type=str,
                        help="a data file for comparing")

    parser.add_argument("-f", "--first", type=str,
                        help="specify the first output file to use")
    parser.add_argument("-s", "--second", type=str,
                        help="specify the second output file to use")

    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="don't print messages")
    parser.add_argument("-o", "--output", type=str, default='lorf2.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")

    parser.add_argument("-L", "--list_genes", action='store_true',
                        help="list common significant genes")
    parser.add_argument("-A", "--findall", action='store_true',
                        help="find all significant concordant genes between given datasets")
    parser.add_argument("-R", "--reversepolarity", action='store_true',
                        help="reverse the polarity of log fold change for second set")



    return parser

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
                bkgd_freq=0.5, label1="group1", label2="group2" ):

    N = 2
    ind = np.arange(2)      # the x locations for the groups
    width = 0.75            # the width of the bars: can also be len(x) sequence

    totals = np.array([float(sum([con_sig1, con_nsig1, ncon_sig1, ncon_nsig1])),
                       float(sum([con_sig2, con_nsig2, ncon_sig2, ncon_nsig2]))])
    verbalise("R", totals)
    ra1 = np.array([con_nsig1, con_nsig2])
    ra2 = np.array([con_sig1, con_sig2])
    ra3 = np.array([ncon_nsig1, ncon_nsig2])
    ra4 = np.array([ncon_sig1, ncon_sig2])


    p1 = plt.bar(ind, ra1/totals, width, color='blue', alpha=0.35)
    p2 = plt.bar(ind, ra2/totals, width, color='blue', alpha=0.8, bottom=(ra1)/totals)
    p3 = plt.bar(ind, ra3/totals, width, color='red', alpha=0.35, bottom=(ra1+ra2)/totals)
    p4 = plt.bar(ind, ra4/totals, width, color='red', alpha=1, bottom=(ra1+ra2+ra3)/totals)
    plt.plot(np.ones(3)*bkgd_freq, color='black', linewidth=2, linestyle='--')
    print "background =", bkgd_freq, np.ones(3)*bkgd_freq
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
                    ha='center', va='bottom')

    autolabel(p1, (ra1)/totals,             ra1/totals, ra1)
    autolabel(p2, (ra1+ra2)/totals,         ra2/totals, ra2)
    autolabel(p3, (ra1+ra2+ra3)/totals,     ra3/totals, ra3)
    autolabel(p4, (ra1+ra2+ra3+ra4)/totals, ra4/totals, ra4)

    plt.show()

def compare_two(args):
    df1 = pd.read_csv(args.first, sep=" ", header=0)
    df2 = pd.read_csv(args.second, sep=" ", header=0)
    df3 = df1.join(df2, lsuffix="1st").dropna()

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

    draw_circles(pos1_u,pos2_u,neg2_u,neg1_u,concord_p,discord_2p,concord_n,discord_1p,label1+" pos", label2+" pos", label1+" neg", label2+" neg")

    print "CONCORDANCE CHARTS"
    #(con_sig1, con_nsig1, ncon_sig1, ncon_nsig1,
    #con_sig2, con_nsig2, ncon_sig2, ncon_nsig2,
    #bkgd) = concordancecounts(df1_sp, df2_sp, df1_sn, df2_sn,
    #                          df1_nsp, df2_nsp, df1_nsn, df2_nsn)

    concordance_sets = concordancecounts(df1_sp, df2_sp, df1_sn, df2_sn,
                              df1_nsp, df2_nsp, df1_nsn, df2_nsn)

    print df3.shape
    draw_graph( *[len(s) for s in concordance_sets[:-1]],
                bkgd_freq=concordance_sets[-1], label1=label1, label2=label2 )


def draw_circles(c1t,c2t,c3t,c4t,o1t,o2t,o3t,o4t,l1t,l2t,l3t,l4t):

    root = tk.Tk()
    root.title('A Circle')

    cw = 600
    ch = 560

    canvas_1 = tk.Canvas(root, width=cw, height=ch, background='white')
    canvas_1.grid(row=0, column=1)

    # specify top left and bottom right coords:
    c1 = (20,40,220,240)
    c2 = (162,40,362,240)
    c3 = (20,182,220,382)
    c4 = (162,182,362,382)

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

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    if args.first and args.second:
        compare_two(args)

    if args.experiments:
        concordant_sig_genesets = []
        for (exp1, exp2) in itertools.product(args.experiments, args.experiments):
            # itertools.product produces an all-by-all comparison of the two lists
            # because this will include the same file compared to itself, we need
            # to remove those cases (which would mess up our comparisons):
            if exp1 == exp2:
                continue

            df1 = pd.read_csv(exp1, sep=" ", header=0)
            df2 = pd.read_csv(exp2, sep=" ", header=0)
            df3 = df1.join(df2, lsuffix="1st").dropna()

            label1 = os.path.basename(exp1)[:5]
            label2 = os.path.basename(exp2)[:5]

            # get genes that are significant in both experiements and their direction:
            df1_sp = set(df3[(df3.padj1st<=0.05) & (df3.log2FoldChange1st>0)].index)
            df2_sp = set(df3[(df3.padj<=0.05) & (df3.log2FoldChange>0)].index)
            df1_sn = set(df3[(df3.padj1st<=0.05) & (df3.log2FoldChange1st<0)].index)
            df2_sn = set(df3[(df3.padj<=0.05) & (df3.log2FoldChange<0)].index)

            df1_nsp = set(df3[(df3.padj1st>0.05) & (df3.log2FoldChange1st>0)].index)
            df2_nsp = set(df3[(df3.padj>0.05) & (df3.log2FoldChange>0)].index)
            df1_nsn = set(df3[(df3.padj1st>0.05) & (df3.log2FoldChange1st<0)].index)
            df2_nsn = set(df3[(df3.padj>0.05) & (df3.log2FoldChange<0)].index)

            print label1, "vs", label2, ":"
            verbalise("G", label1, "pos", len(df1_sp))
            verbalise("G", label1, "neg", len(df1_sn))
            verbalise("G", label1, "both", len(df1_sp | df1_sn))
            verbalise("C", label2, "pos", len(df2_sp))
            verbalise("C", label2, "neg", len(df2_sn))
            verbalise("C", label2, "both", len(df2_sp | df2_sn))

            concordance_sets = concordancecounts(df1_sp, df2_sp, df1_sn, df2_sn,
                              df1_nsp, df2_nsp, df1_nsn, df2_nsn)

            concordant_sig_genesets.append(concordance_sets[0])

            verbalise("Y", "concordant DEGs = ", len(concordance_sets[0]))

        common_to_all = set.intersection(*concordant_sig_genesets)
        verbalise("R", "There are %d genes common to all datasets" % len(common_to_all))

        if args.list_genes:
            ncbi = bm.ncbi_dic()
            print "\n".join([bm.show_name(loc, ncbi) for loc in common_to_all])
