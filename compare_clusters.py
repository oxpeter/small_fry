#!/usr/bin/python

"""
module used to compare gene clusters between forg to stat and stat to forg
transitions.
"""
import os
import sys
import re
import argparse
import itertools
import string
import random

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

from genomepy import config

############################################
def define_arguments():
    parser = argparse.ArgumentParser(description=
            "Compares the gene membership of multiple lists within two directories")

    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='genematch.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")

    # data file options:
    parser.add_argument("input", type=str, nargs='+',
                        help="directories containing cluster lists to compare")
    parser.add_argument("-c", "--column", type=int, default=0,
                        help="column in which gene names are found (default = 0)")
    parser.add_argument("-t", "--filetype", type=str, default='list|txt',
                        help="""specify the file filter to apply to each directory
                        (default = list|txt)""")
    parser.add_argument("-f", "--filter", type=str, default=None,
                        help="""[FILE,COLUMN] Create a filter list from column COLUMN
                        in file FILE. Comparisons will then only look at genes in each
                        cluster that are also found in the file specified (useful for
                        looking at only significantly differentially expressed genes)""")
    parser.add_argument("-w", "--weight", type=int, default=2,
                        help="""Threshold number of items to be shared by two clusters
                        for a solid line to be drawn between them in the network. If there
                        are less than this number, a dotted line will be drawn instead""") 

    parser.add_argument("-n", "--network", type=int, default=0,
                        help="""Calculate and display a weighted network graph. 
                        1: circular, 2: random, 3: spectral, 4: spring, 5: shell,
                        6:graphviz""")

    return parser

def jaccard(s0, s1):
    "returns the Jaccard index of two sets"
    denom = len((s0 | s1))
    if denom == 0:
        return -1.
    else:
        return 1. * len((s0 & s1)) / denom


def update_dic(dic, name, score, partner, jidx, reciprocal=True):
    assert isinstance(score, int), "score is not an integer: %r" % score
    assert isinstance(jidx, float), "jaccard idx is not a float: %r" % jidx
    if name in dic:
        if score > dic[name][1]:
            dic[name] = (partner, score, jidx)
    else:
        dic[name] = (partner, score, jidx)
    if reciprocal:
        if partner in dic:
            if score > dic[partner][1]:
                dic[partner] = (name, score, jidx)
        else:
            dic[partner] = (name, score, jidx)

def filter_list(list, filterlist, filter=True):
    if filter:
        return [ i for i in list if i in filterlist ]
    else:
        return list

def clean_filename(filename):
    fsearch = re.search("([FS]2[FS]).cluster_([0-9]+)",filename)
    if fsearch:
        return fsearch.group(1) + fsearch.group(2)
    else:
        return "?"

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)
    outfile = logfile[:-3] + "out"

    if args.filter:
        filterset = config.make_a_list(args.filter.split(',')[0], int(args.filter.split(',')[1]))
    else:
        filterset = []

    # collect lists of files to compare:
    file_lists = {}
    verbalise("B", "Collecting lists of files using %s filter" % args.filetype)
    for i, dir in enumerate(args.input):
        file_lists[i] = [ os.path.join(dir,f) for f in os.listdir(dir) if re.search(args.filetype, f) ]
        verbalise("Y", "%d files found for dir %s" % (len(file_lists[i]), dir))
    # compare each set and save best matches to dictionary best_match
    best_match = {}
    network_weights = {}
    collect_clusters = []
    for pair in itertools.permutations(file_lists.keys(), 2):
        for cpair in itertools.product(file_lists[pair[0]], file_lists[pair[1]]):
            s0n = cpair[0]
            s1n = cpair[1]
            s0  = set(filter_list(config.make_a_list(cpair[0], args.column),filterset, args.filter))
            s1  = set(filter_list(config.make_a_list(cpair[1], args.column),filterset, args.filter))
            score = len(s0 & s1)
            jidx = jaccard(s1, s0)
            update_dic(best_match, s0n, score, s1n, jidx)
            if score != 0:
                if s0n in network_weights:
                    network_weights[s0n][s1n] = score
                else:
                    network_weights[s0n] = {s1n:score}
                if s1n in network_weights:
                    network_weights[s1n][s0n] = score
                else:
                    network_weights[s1n] = {s0n:score}

    # write best matches to file:
    handle = open( outfile, 'w')
    for file in best_match:
        handle.write( "%s %s %d %.2f\n" % (os.path.basename(file),
                                os.path.basename(best_match[file][0]),
                                best_match[file][1], best_match[file][2]))
    handle.close()

    verbalise("Y", "FILE1\tFILE2\t# matches\tJaccard Index")
    os.system("sort -k 4,4n " + outfile)

    if args.network:
        df = pd.DataFrame(network_weights).fillna(0)
        cluster_array = df.values
        dt = [('len', float)]
        cluster_array = cluster_array.view(dt)
        
        # create networkx object:
        G = nx.from_numpy_matrix(cluster_array)
        G = nx.relabel_nodes(G, dict(zip(range(len(G.nodes())),[ clean_filename(f) for f in df.columns])))
        
        # add weights to each edge for later processing:
        e = [ (n1, n2, G[n1][n2]['len']) for n1 in G for n2 in G[n1]   ]
        G.add_weighted_edges_from(e)

        # define the type of graph to be drawn:        
        network_types = {1: nx.circular_layout, 
                    2: nx.random_layout, 
                    3: nx.spectral_layout, 
                    4: nx.spring_layout, 
                    5: nx.shell_layout, 
                    6: nx.graphviz_layout}
        net_type = network_types[args.network]
        
        # check for help in parsing the network objects:
        if 'S2F77' in G:
            verbalise( "G['S2F77'] --> %s" % (G['S2F77']))
        
        # split into all sub-networks and draw each one:
        pos=net_type(G)
        C=nx.connected_component_subgraphs(G)

        for g in C:
            # report size of sub-network
            verbalise("G", "%d clusters in sub-network" % (len(g)) )
            
            # define which edges are drawn bold:
            elarge=[(u,v) for (u,v,d) in g.edges(data=True) if d['weight'] >=args.weight]
            esmall=[(u,v) for (u,v,d) in g.edges(data=True) if d['weight'] < args.weight]
            verbalise("M", "%d cluster pairs with more than %d shared members" % (len(elarge),args.weight)) 
            verbalise("G", "\n".join([str(t) for t in g.edges(data=True) if d['weight'] >=args.weight]))
            
            # draw edges:
            nx.draw_networkx_edges(g,pos,edgelist=elarge,
                    width=2)
            nx.draw_networkx_edges(g,pos,edgelist=esmall,
                    width=1, alpha=0.3, edge_color='blue', style='dashed')
        
            # draw sub-network:
            nx.draw(g,
                 pos,
                 node_size=40,
                 node_color=[ {'F':1, 'S':2}[n[0]] for n in g ],
                 vmin=0.0,
                 vmax=2.0,
                 width=0,
                 with_labels=True
                 )
        
        plt.show()
        



