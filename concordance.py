#!/usr/bin/env python
"""
This program is a simulation of gene expression studies performed in four different ant species,
in order to determine the likelihood of achieving given congruence between studies by chance.
The starting point is the Cerapachys biroi, Acromyrmex echinatior, Solenopsis invicta and
Dinoponera quadriceps queen/worker RNA-seq experiements, which show only 5 genes concordantly
differentially expressed between all four species.

The program will take 7996 'genes' that are equivalent to the 7996 unambiguously orthologous
genes in the four species, and randomly assign equivalent numbers of them to be significant
DEGs to each species, then plot the degree of overlap.

Later I hope to also incoporate some sort of correlation of gene expression as a function of
the distance between species.
"""

import argparse
import os
import tempfile
import itertools
import random
import re
import datetime

import matplotlib.pyplot as plt
import numpy as np

from ortholotree import config

############################################################################

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "This program performs subset analysis of randomly produced DEG lists")
    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="Print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='concordance',
                        help="Specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="Specify the directory to save results to")
    parser.add_argument("-D", "--display_on", action='store_true',default=False,
                        help="Display graph results (eg for p value calculation)")
    parser.add_argument("-A", "--save_all", action='store_true',default=False,
                        help="Save all graph results from each iteration [default = False]")

    # analysis options:
    parser.add_argument("-i", "--iterations", type=int, default=100,
                        help="Specify the number of iterations to run [default = 100]")
    parser.add_argument("-t", "--total_orthos", type=int, default=5008,
                        help="""Specify the total number of orthologs from which to select
                        the differentially expressed genes [default = 5008] """)
    parser.add_argument("-s", "--degs", type=str, default='2156,2938,1795,306',
                        help="""Specify a comma-separated list of numbers representing
                        the number of DEGs present in each sample to be compared.
                        [default = '2156,2938,1795,306']""")
    parser.add_argument("--significance", type=str,
                        help="""Calculate the significance of specified numbers from
                         the distribution of all-species overlap""")
    parser.add_argument('-c', "--concordance", type=float,
                        help="""a number between 0 and 1 that indicates the probability of
                        a gene being concordant between two species. Default = 1, meaning
                        all significant genes are concordant.""")


    return parser

def assign(concordance=1):
    """
    concordance is the probability that two genes have the same polarity.
    the probability of assigning a direction to a gene, such that this concordance
    probability is met, is determined from the quadratic equation.

    Pr(concordant) = p^2 + q^2   where p is the probability of being assigned +ve
                                and q is the probabliity of being assigned -ve
    (p + q) = 1

    therefore, Pr(concordant) = 1 - 2p + 2p^2

    or,        2p^2 - 2p + 1 - Pr(conc) = 0
    therefore                         p = (2 - sqrt(4 - 8 * (1-Pr)) ) / 4
    """
    # calculate p, the probability that a gene is positive:
    p = (2 - np.sqrt(4 - 8 * (1-concordance)) ) / 4

    # calculate a random number, then see if it is greater than or less than 1000p
    marker = random.randint(0,999)
    if marker >= 1000 * p:
        return '+'
    else:
        return '-'

def make_genes(numgenes=7996, concordance=1):
    genelist = [ "%s%d" % (assign(concordance), i) for i in range(numgenes) ]
    return genelist

def select_degs(genelist, numgenes=1000):
    return set(random.sample(genelist, numgenes))

def find_overlaps(*degsets):
    """
    This will work through all the possible sets of overlap between 2, 3, 4... sets.
    will return { 1:[x1,x2,x3,x4], 2:[y1, y2,...], 3:[z1,z3,z3,z4] ... } where the key is
    the number of overlapping sets, and the value is a list of the size of all
    possible combinations of that number of sets.
    """
    common_degs = {}
    # determine number of common genes in 1 to n-sized venn diagrams
    for compsize in range(1,len(degsets) + 1):
        for i,deg_combo in enumerate(itertools.combinations(degsets,compsize)):
            # deg_combo is a list of size compsize
            common_degs[(i, len(deg_combo))] = set.intersection(*deg_combo)


    # determine how many unique items in each part of an n-way venn diagram:
    nos = {} # non-overlapping sets for each cateogory of overlap (common to 2,3,4 etc sets)
    for degcombo in common_degs:
        if degcombo[1] not in nos:
            nos[degcombo[1]] = []
        unique_set = reduce(lambda x,y: x - y,
                [common_degs[degcombo]] + [ common_degs[k] for k in common_degs.keys() if k[1] > degcombo[1]  ] )
        nos[degcombo[1]].append(len(unique_set))

    return nos

def add_wobble(X):
    "Adds noise to enable easier viewing of multiple points"
    newX = [ x - 1.0/8 + random.random() / 4 for x in X ]
    return newX

def graph_nos(nos, outfile, display=False, save=False, title=None):
    # create boxplot with overlapping individual data points
    points = [ (x,y) for x in nos for y in nos[x] ]
    annocoords = [ (x,np.mean(nos[x])) for x in nos ]
    annomeans  = [ "%.1f\n( +/- %.1f )" % (np.mean(nos[x]),
                                np.std(nos[x])/np.sqrt(np.mean(nos[x]))) for x in nos ]
    X, Y = zip(*points)
    if len(X) > 250:
        X = add_wobble(X)

    plt.boxplot(nos.values())
    plt.plot(X, Y, 'co', alpha=0.5)

    # annotate graph
    if not title:
        plt.title("Number of genes common to n overlapping sets")
    else:
        plt.title(title)
    plt.xlabel('number of overlapping sets')
    plt.ylabel('number of common genes')

    for idx in range(len(annocoords)):
        plt.annotate(annomeans[idx], xy=annocoords[idx], xycoords='data',
                    xytext=(35, 0), textcoords='offset points',
                    horizontalalignment='left', verticalalignment='bottom',
                    )
    # polish and publish
    plt.tight_layout()
    if save:
        plt.savefig(outfile[:-3] + "png", format='png')
    if display:
        plt.show()
    else:
        plt.close()

def graph_final(all_values, outfile, display=False, save=False, title=None):
    # create boxplot with overlapping individual data points

    annocoords = (1.05,np.mean(all_values))
    annomeans  = "%.1f\n( +/- %.1f )" % (np.mean(all_values), np.std(all_values))

    plt.figure(figsize=(5,4), facecolor='blue')
    plt.boxplot(all_values)
    plt.scatter(1, np.mean(all_values), marker='*')
    # annotate graph
    if not title:
        plt.title("Number of genes common to all overlapping sets")
    else:
        plt.title(title)
    plt.ylabel('number of concordant DEGs')

    plt.annotate(annomeans, xy=annocoords, xycoords='data',
                xytext=(35, 0), textcoords='offset points',
                horizontalalignment='left', verticalalignment='bottom',
                )
    # polish and publish
    plt.tight_layout()
    if save:
        plt.savefig(outfile[:-3] + "png", format='png')
    if display:
        plt.show()
    else:
        plt.close()

def report_all(outfile, all_values, significance_values):
    handle = open(outfile, 'w')
    for k in all_values.keys():
        for s in [ int(r) for r in significance_values]:
            pval = 1.0 * sum( 1 for x in all_values[k] if x >= s ) / len(all_values[k])
            pval_str = "P-value for %d genes common to %d species is %.3f" % (s, k, pval)
            verbalise("G", pval_str)
            handle.write(pval_str + "\n")

        distrib = "%.2f +/- %.2f [%d - %d]" % (np.mean(all_values[k]),
                                                np.std(all_values[k]),
                                                min(all_values[k]),
                                                max(all_values[k]))
        verbalise("Y", distrib)
        handle.write(distrib + "\n")
    handle.close()

def report_final(outfile, all_values, significance_values, k):
    handle = open(outfile, 'w')
    for s in [ int(r) for r in significance_values]:
        pval = 1.0 * sum( 1 for x in all_values if x >= s) / len(all_values)
        pval_str = "P-value for %d genes common to %d species is %.3f" % (s, k, pval)
        verbalise("G", pval_str)
        handle.write(pval_str + "\n")

        distrib = "%.2f +/- %.2f [%d - %d]" % (np.mean(all_values),
                                                np.std(all_values),
                                                min(all_values),
                                                max(all_values))
        verbalise("Y", distrib)
        handle.write(distrib + "\n")
    handle.close()

############################################################################

if __name__ == '__main__':
    dbpaths = config.import_paths()

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    temp_dir = tempfile.mkdtemp()

    means = {}
    everything = {}
    final_count = []
    for iter in range(args.iterations):
        degs = {}
        # for each sample, create a genome and select the indicated number of
        # DEGs from it, adding the list to the dictionary degs[sample] = [DEGs]
        for i,numdegs in enumerate([ int(i) for i in args.degs.split(',') ]):
            genome = make_genes(numgenes=args.total_orthos,
                                    concordance=args.concordance)

            degs[(i,numdegs)] = select_degs(
                                    genome,
                                    numdegs)

        # perform breakdown analysis:
        venn_values = find_overlaps(*degs.values())
        for k in venn_values:
            if k not in means:
                means[k] = []
            if k not in everything:
                everything[k] = []
            means[k].append(np.mean(venn_values[k]))
            for v in venn_values[k]:
                everything[k].append(v)

        if args.save_all:
            graph_nos(venn_values,
                    logfile[:-4] + "_iter%d.png" % iter,
                    display=args.save_all and args.display_on,
                    save=args.save_all)

        # get final number of concordant DEGs:
        common_degs = set.intersection(*[set(l) for l in degs.values()])
        final_count.append(len(common_degs))


    verbalise("M", "stats for last genome generated:")
    verbalise("Y",
    "Genome size: %d\nNum DEGs: %d\nRatio + to -: %.3f\nFirst 10 genes:\n%s\n\n" % (len(genome),
                            numdegs,
                            1.*sum(1 for i in genome if i[0]=='+')/len(genome),
                            " ".join(genome[:10]),))


    if args.significance:
        report_all(logfile[:-3] + "subset_pvalues.out",
                    everything,
                    args.significance.split(','))
        report_final(logfile[:-3] + "pvalues.out",
                    final_count,
                    args.significance.split(','),
                    k=len(args.degs.split(',')))

    # generate graphs:
    graph_final(final_count,
                    logfile[:-3] + "all_concordant_iterations.png",
                    display=args.display_on,
                    save=True,
                    title = "# concordant DEGs, %d spp\n(%d iterations, P(conc)=%.2f)" % (len(args.degs.split(',')),args.iterations, args.concordance))


    graph_nos(means,
                    logfile[:-3] + "mean.png",
                    display=args.display_on,
                    save=True,
                    title = "mean (+/- SEM) number of concordant random DEGs \nfrom %d iterations, Pr(concordant)=%.2f" % (args.iterations, args.concordance))

    graph_nos(everything,
                    logfile[:-3] + "all.png",
                    display=args.display_on,
                    save=True,
                    title = "all values of %d iterations" % args.iterations)

    fake_graph = dict.fromkeys(range(max(everything.keys())),[1])
    fake_graph[max(everything.keys())] = everything[max(everything.keys())]
    graph_nos(fake_graph,
                    logfile[:-3] + "all_max_overlap.png",
                    display=args.display_on,
                    save=True,
                    title = "mean values of %d iterations" % args.iterations)


