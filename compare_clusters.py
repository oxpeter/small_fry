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
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist

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
    parser.add_argument("input", metavar='<CLUSTER_DIR>', type=str, nargs='+',
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
    parser.add_argument("--second_degree", action='store_true',default=False,
                        help="analyse which clusters are separated by 2 degrees of separation")
    parser.add_argument("--jidx", action='store_true',default=False,
                        help="""Use jaccard index for drawing network instead of item
                            number. If set, make sure the weights specified are values
                            between 0 and 1.""")

    # viewing options:
    parser.add_argument("--display_off", action='store_true',default=False,
                        help="Do not display plots (images will still be saved to pdf)")

    parser.add_argument("-w", "--minweight", type=float, default=2,
                        help="""The minimum number of items to be shared by two clusters
                        for a solid blue line to be drawn between them in the network. If there
                        are less than this number, a blue dotted line will be drawn instead
                        (default = 2)""")
    parser.add_argument("-W", "--maxweight", type=float, default=10,
                        help="""The minimum number of items to be shared by two clusters
                        for a thick solid black line to be drawn between them in the
                        network. If there are less than this number, a solid blue line will
                        be drawn instead (default = 10)""")
    parser.add_argument("-n", "--network", type=int, default=0,
                        help="""Calculate and display a weighted network graph.
                        1: circular, 2: random, 3: spectral, 4: spring, 5: shell,
                        6:pygraphviz""")
    parser.add_argument('-D', "--dimensions", type=int, default=2,
                        help="chose between 2 or 3 dimensions")

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
    fsearch = re.search("([FS]2[FS]).cluster_([0-9]+)",filename) # for the FS24 cluster names
    csearch = re.search("cluster(\w+)",filename) # cluster search to find cluster id
    gsearch = re.search("(\d+)", os.path.basename(filename)) # very generic search to pull a unique id
    if fsearch:
        return fsearch.group(1)[0] + fsearch.group(1)[-1] + ' ' + fsearch.group(2)
    elif csearch:
        return csearch.group(1)
    elif gsearch:
        return os.path.basename(filename)[0:5] + "_" + gsearch.group(1)
    else:
        return os.path.basename(filename)[0:5]

def collect_weights(file_lists, store_jidx=False):
    """
    INPUT:
    file_lists = { 0: [(name1,list1),(name2, list2)], 1:[(name3,list3),(name4, list4)]  }
    store_jidx -> if False, then score will be appended to dictionary. if True, then jidx.
    """
    # compare each set and save best matches to dictionary best_match
    best_match = {}
    network_weights = {}    # to collect the number of items shared between a pair of clusters
    cluster_sizes = {}      # to collect the total number of items in each cluster.

    if len(file_lists.keys()) > 1:
        dir_iter = itertools.permutations(file_lists.keys(), 2)
    elif len(file_lists.keys()) == 1:
        dir_iter = ((0, 0),)
    else:
        verbalise("R", "No list supplied!! Exiting")
        exit()

    for dir_pair in dir_iter:
        for cpair in itertools.product(file_lists[dir_pair[0]], file_lists[dir_pair[1]]):
            s0n = cpair[0][0]
            s1n = cpair[1][0]
            if s0n == s1n:
                continue
            s0  = set(cpair[0][1])
            s1  = set(cpair[1][1])
            score = len(s0 & s1)
            jidx = jaccard(s1, s0)
            update_dic(best_match, s0n, score, s1n, jidx)
            cluster_sizes[s0n] = len(s0)
            cluster_sizes[s1n] = len(s1)

            if store_jidx:
                val = jidx
            else:
                val = score

            if score != 0:
                if s0n in network_weights:
                    network_weights[s0n][s1n] = val
                else:
                    network_weights[s0n] = {s1n:val}
                if s1n in network_weights:
                    network_weights[s1n][s0n] = val
                else:
                    network_weights[s1n] = {s0n:val}

    return best_match, network_weights, cluster_sizes

def draw_network(network_weights, maxweight=20, minweight=10, dims=2, report=None, display=True):
    """
    Given a dictionary of weights between nodes, draws the network structure.

    input:
        node1[node2] = weight

    """
    df = pd.DataFrame(network_weights).fillna(0)
    cluster_array = df.values
    dt = [('len', float)]
    cluster_array = cluster_array.view(dt)

    # create networkx object:
    G = nx.from_numpy_matrix(cluster_array)
    relabel = dict(zip(range(len(G.nodes())),[ clean_filename(f) for f in df.columns]))
    if dims == 2:
        G = nx.relabel_nodes(G, relabel)

    # create dictionary to convert names back to positional labels:
    #backlabels = dict(zip([ clean_filename(f) for f in df.columns], range(len(G.nodes()))))
    #print backlabels.items()[:5]

    # add weights to each edge for later processing:
    e = [ (n1, n2, G[n1][n2]['len']) for n1 in G for n2 in G[n1]   ]
    G.add_weighted_edges_from(e)

    # define the type of graph to be drawn:
    network_types = {1: nx.circular_layout,
                2: nx.random_layout,
                3: nx.spectral_layout,
                4: nx.spring_layout,
                5: nx.shell_layout,
                6: nx.pygraphviz_layout}
    net_type = network_types[args.network]

    # check for help in parsing the network objects:
    if 'S2F99' in G:
        verbalise( "G['S2F99'] --> %s" % (G['S2F99']))

    # split into all sub-networks and draw each one:
    if dims == 3:
        from mayavi import mlab
        pos=net_type(G, dim=dims, k=0.15)
    else:
        pos=net_type(G)
    C = nx.connected_component_subgraphs(G)

    # initialise pdf for saving all plots:
    if report:
        pp = PdfPages( report[:-3] + 'pdf' )
    for g in C:
        # report size of sub-network
        verbalise("Y", "%d clusters in sub-network" % (len(g)) )

        # define which edges are drawn bold:
        rlarge =  [(u,v,d) for (u,v,d) in g.edges(data=True) if d['weight'] >= maxweight]
        rmedium =[(u,v,d) for (u,v,d) in g.edges(data=True) if maxweight > d['weight'] >= minweight]
        rsmall =  [(u,v,d) for (u,v,d) in g.edges(data=True) if d['weight'] < minweight]

        elarge =  [ (u,v) for (u,v,d) in rlarge  ]
        emedium = [ (u,v) for (u,v,d) in rmedium ]
        esmall =  [ (u,v) for (u,v,d) in rsmall  ]

        rlarge.sort(key=lambda x: x[2]['weight'])
        rmedium.sort(key=lambda x: x[2]['weight'])
        rsmall.sort(key=lambda x: x[2]['weight'])

        # report number of clusters with each weight
        verbalise("M", "%d cluster pairs with %d or more shared members" % (len(elarge),maxweight))
        verbalise("G", "\n".join([ "%-6r %-6r %s" % (t[0], t[1], t[2]['weight']) for t in rlarge]))

        verbalise("M", "%d cluster pairs with less than %d and %d or more shared members" % (len(emedium),maxweight, minweight))
        verbalise("G", "\n".join(["%-6r %-6r %s" % (t[0], t[1], t[2]['weight']) for t in rmedium][-3:]))

        verbalise("M", "%d cluster pairs with less than %d shared members" % (len(esmall),minweight))
        verbalise("G", "\n".join(["%-6r %-6r %s" % (t[0], t[1], t[2]['weight']) for t in rsmall][-3:]))
        verbalise("G","")

        if report:
            handle = open(report, 'a')

            handle.write("%d clusters in sub-network\n" % (len(g)))
            handle.write("%d cluster pairs with %d or more shared members\n" % (len(elarge),maxweight))
            handle.write("%d cluster pairs with less than %d and %d or more shared members\n" % (len(emedium),maxweight, minweight))
            handle.write("%d cluster pairs with less than %d shared members\n" % (len(esmall),minweight))
            handle.write("\n".join(["%-6r %-6r %s" % (t[0], t[1], t[2]['weight']) for t in rlarge]) + '\n')
            handle.write("\n".join(["%-6r %-6r %s" % (t[0], t[1], t[2]['weight']) for t in rmedium]) + '\n')
            handle.write("\n".join(["%-6r %-6r %s" % (t[0], t[1], t[2]['weight']) for t in rsmall]) + '\n\n')
            handle.close()

        if dims == 2:
            # draw edges (partitioned by each edge type):
            # large:
            nx.draw_networkx_edges(g,pos,edgelist=elarge,
                    width=2, edge_color='purple')
            # medium:
            nx.draw_networkx_edges(g,pos,edgelist=emedium,
                    width=1, alpha=0.6, edge_color='blue')
            # small:
            nx.draw_networkx_edges(g,pos,edgelist=esmall,
                    width=1, alpha=0.3, edge_color='blue', style='dashed')

            # draw sub-network:
            nx.draw(g,
                 pos,
                 node_size=40,
                 node_color=[ {'F':'g', 'S':'b', '_':'r'}[n[0]] for n in g ],
                 vmin=0.0,
                 vmax=2.0,
                 width=0,
                 with_labels=True
                 )
            if report:
                plt.savefig(pp, format='pdf')
            if display:
                plt.show()
            else:
                plt.close()

        if dims == 3:

            xyz=np.array([pos[n] for n in g ])

            mlab.figure(1, bgcolor=(0, 0, 0))
            mlab.clf()
            for shape, edges, colour in [('sphere', elarge,  (0.5,0.9,0.2)),
                                         ('sphere', emedium, (0.2,0.5,0.9)),
                                         ('sphere', esmall,  (0.9,0.2,0.5)),
                                        ]:

                pts = mlab.points3d(xyz[:,0], xyz[:,1], xyz[:,2],
                    [ {'F':(4), 'S':(6), '_':(8)}[relabel[n][0]] for n in g ],
                    colormap = 'copper',
                    scale_factor=0.1,
                    scale_mode='none',
                    resolution=20,
                    mode=shape,
                    vmax=10,
                    vmin=0)

                edge_array =  np.array([ list(t) for t in edges ])
                pts.mlab_source.dataset.lines = edge_array
                tube = mlab.pipeline.tube(pts, tube_radius=0.01)
                mlab.pipeline.surface(tube, color=colour, vmin=0, vmax=10)
            if display:
                mlab.show()

    if report:
        pp.close()
    return rlarge, rmedium, rsmall

def plot_heatmap(weights, outfile=None, samplesx=None, samplesy=None):
    """
    INPUT:
    weights -> a dictionary of values for each XY pair: weights[X][Y] = w
    outfile -> file to save heatmap image to
    samplex -> list of sample names for X axis. If None, all samples will be used.
    sampley -> list of sample names for Y axis. If None, all samples will be used
    """
    df = pd.DataFrame( weights )

    if not samplesx:
        samplesx = samplesy
    elif not samplesy:
        samplesy = samplesx

    if samplesx:
        selectx = [ name for name in samplesx if name in df.columns ]
        selecty = [ name for name in samplesy if name in list(df.index.values) ]
        namesx = [ clean_filename(name) for name in samplesx if name in df.columns ]
        namesy = [ clean_filename(name) for name in samplesy if name in list(df.index.values) ]

        df = df[selectx].loc[selecty]
        X = np.nan_to_num(df.values)
    else:
        X = np.nan_to_num(df.values)
        namesx = [ clean_filename(name) for name in df.columns ]
        namesy = [ clean_filename(name) for name in df.columns ]

    #verbalise("B", df.columns)
    #verbalise("C", namesx)
    #verbalise("Y", namesy)

    # plot top dendrogram
    fig = plt.figure(figsize=(8, 8))
    axx = fig.add_axes([0.22,0.76,0.6,0.2], frame_on=True)
    dx = dist.pdist(X.T)
    Dx = dist.squareform(dx)
    Yx = sch.linkage(Dx, method='complete', metric='euclidean')
    #np.clip(Yx[:,2], 0, 100000, Yx[:,2]) # prevents errors from negative floating point near zero
    Zx = sch.dendrogram(Yx)

    # plot left dendrogram
    axy = fig.add_axes([0.01,0.15,0.2,0.6], frame_on=True)
    dy = dist.pdist(X)
    Dy = dist.squareform(dy)  # full matrix
    Yy = sch.linkage(Dy, method='complete', metric='euclidean')
    Zy = sch.dendrogram(Yy, orientation='right')

    # remove ticks from dendrograms
    axx.set_xticks([])
    axx.set_yticks([])
    axy.set_xticks([])
    axy.set_yticks([])

    # reorder matrices
    indexx = Zx['leaves']
    indexy = Zy['leaves']
    Dy = Dy[indexy,:]
    Dy = Dy[:,indexx]
    X = X[indexy,:]
    X = X[:,indexx]
    newnamesx = [ namesx[i] for i in indexx ]
    newnamesy = [ namesy[i] for i in indexy ]

    # display distance matrix
    axmatrix = fig.add_axes([0.22,0.15,0.6,0.6])
    im = axmatrix.matshow(X, aspect='auto', origin='lower')

    # set position and naming of axes
    axmatrix.set_xticks(range(len(newnamesx)))
    axmatrix.set_yticks(range(len(newnamesy)))
    axmatrix.set_xticklabels(newnamesx, rotation=90)
    axmatrix.set_yticklabels(newnamesy)
    axmatrix.xaxis.tick_bottom()
    axmatrix.yaxis.tick_right()
    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90)

    if outfile:
        plt.savefig(outfile)
    plt.show()

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)
    outfile = logfile[:-3] + "out"

    # check that weights and scoring method are compatible (and meaningful):
    if args.jidx and ( args.minweight > 1 or args.maxweight > 1 ):
        verbalise("R", """You have selected the jaccard index for measuring distance
        between  samples, but have selected values greater than 1 for your minimum and
        maximum weights to display in the network graph.""")
        exit()

    if args.filter:
        filterset = config.make_a_list(args.filter.split(',')[0], int(args.filter.split(',')[1]))
    else:
        filterset = []

    # collect lists of genes to compare:
    file_lists = {}
    verbalise("B", "Collecting lists of files using %s filter" % args.filetype)
    for i, dir in enumerate(args.input):
        file_lists[i] = [ (os.path.join(dir,f),
                            filter_list(config.make_a_list(os.path.join(dir,f),args.column),
                                        filterset,
                                        args.filter)) for f in os.listdir(dir) if re.search(args.filetype, f) ]
        verbalise("Y", "%d files found for dir %s" % (len(file_lists[i]), dir))


    # compare each set and save best matches to dictionary best_match
    best_match, network_weights, cluster_sizes = collect_weights(file_lists, store_jidx=args.jidx)


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
        edge_weights = draw_network(network_weights,
                                    minweight=args.minweight,
                                    maxweight=args.maxweight,
                                    dims=args.dimensions,
                                    report=outfile,
                                    display=(not args.display_off),
                                    )


    if not args.display_off:
        # determine the maximum bin size to fix both histograms to same x-axis range:
        max_x = max(cluster_sizes.values())

        # calculate histograms (separate y-axes)
        fig, ax1 = plt.subplots()
        ax1.hist([ network_weights[p1][p2] for p1 in network_weights for p2 in network_weights[p1]],
                    bins=26,
                    color='red',
                    alpha=0.6,
                    range=[0,max_x])
        ax1.set_xlabel('number of genes shared/\nnumber of genes in clusters', color='black')
        ax1.set_ylabel('number of edges', color='r')

        ax2 = ax1.twinx()
        ax2.hist( cluster_sizes.values(),
                    bins = 26,
                    color='green',
                    alpha=0.6,
                    range=[0,max_x])
        ax2.set_ylabel('number of clusters', color='g')

        plt.show()
    else:
        plt.close()

    # get names of samples from each list of clusters:
    x_list = [ l[0] for l in file_lists[1] if l[0] in network_weights]
    y_list = [ l[0] for l in file_lists[0] if l[0] in network_weights]

    plot_heatmap(network_weights,
                    logfile[:-3] + 'heatmap.pdf',
                    x_list,
                    y_list)


    if args.second_degree:
        # calculate closest clusters from same folder (based on 2 degrees of separation):
        independent_sets = []
        for n1 in network_weights:
            independent_sets.append((n1,
                [ n2 for n2 in network_weights[n1] if network_weights[n1][n2] >0 and n1 != n2 ]))

        comp_dic = {0:independent_sets}
        d2best_match, d2network_weights, d2cluster_sizes = collect_weights(comp_dic, store_jidx=args.jidx)
        if args.network:
            edge_weights = draw_network(d2network_weights,
                                        minweight=args.minweight,
                                        maxweight=args.maxweight,
                                        dims=args.dimensions,
                                        report=logfile[:-3] + "2nd_degree.out",
                                        display=(not args.display_off),
                                        )

        plot_heatmap(d2network_weights,
                    logfile[:-3] + 'heatmap1.pdf',
                    x_list)

        plot_heatmap(d2network_weights,
                    logfile[:-3] + 'heatmap2.pdf',
                    y_list)

