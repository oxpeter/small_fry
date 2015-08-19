#! /usr/bin/env python
"""
A wrapper to identify, extract, align and phylogenise a target gene. The primary purpose
of this module is for a user to confirm orthology of their gene, not to understand in
detail the phylogenetic relationships between all closely related proteins.
"""

import argparse
import linecache
from math import floor
import os
import re
import sys
import tempfile

from Bio import AlignIO
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np

from genomepy import config
from genomepy import genematch as gm

############################################################################

dbpaths = config.import_paths()

############################################################################

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "A module to perform a variety of gene term related analyses")
    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='genematch.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-D", "--display_on", action='store_true',default=False,
                        help="display graph results (eg for p value calculation)")


    # data file options:
    parser.add_argument("-f", "--fasta", type=str,
                        help="Fasta file of gene to analyse (must be amino acids)")
    parser.add_argument("-g", "--gene", type=str,
                        help="""Gene or transcript name. This name will be checked
                        against the peptide database to see if there is a match.""")
    parser.add_argument("-s", "--species", type=str, default="",
                        help="""Four letter abbreviation of species, to speed up
                        searching for genes. E.g Cerapachys biroi = Cbir""")
    parser.add_argument("-n", "--name_conversion", type=str,
                        help="""Name conversion file for RAxML only analysis. Format is
                        the same as that outputted by a standard run. If rerunning an
                        analysis, it is therefore possible to supply the mafft.out
                        file and the name_conversion.txt file to create a phylogeny with
                        the original peptide names.""")

    # analysis options:
    parser.add_argument("-B", "--buildhmmer", action='store_true',
                        help="""Use sequences supplied in fasta file or gene to build a
                        hmmer model and use hmmer to extract sequence matches (more
                        sensitive than using the default phmmer setting, but requires
                        knowledge of the genes to input). This may be useful after a
                        preliminary run, using a resulting orthologous cluster as the
                        new input.""")
    parser.add_argument("-c", "--mincollect", type=int, default=2,
                        help="""minimum number of homologs to collect from each species,
                        regardless of how poor the similarity""")
    parser.add_argument("-t", "--scorethresh", type=float, default=0.5,
                        help="""minimum fraction of the score of the best matching
                        gene for a species for including additional matches""")
    parser.add_argument("-p", "--threads", type=int, default=2,
                        help="number of threads to use for RAxML calculations")
    parser.add_argument("-b", "--bootstrap", type=int,
                        help="""Perform rapid bootstrap analysis. Requires a random
                         integer for the bootstrap seed.""")
    parser.add_argument("-a", "--globalthresh", type=float, default=0.2,
                        help="""Identify the top proportion of all results globally.""")
    parser.add_argument("-R", "--raxml_only", action='store_true',
                        help="""Construct phylogeny using RAxML using fasta alignment
                        provided in fasta file """)
    parser.add_argument("-x", "--exclude_genes", type=str,
                        help="""A comma-separated list of genes to exclude from the
                        alignment and phylogeny""")
    parser.add_argument("-e", "--exclude_species", type=str,
                        help="""A comma-separated list of species to exclude from the
                        alignment and phylogeny. Use the four-letter abbreviation.""")
    parser.add_argument("-l", "--maxlength", type=int,
                        help="""If provided, will remove all genes longer than this
                        size. Useful for removing concatenated genes that otherwise
                        are orthologous. It is not recommended that you use this
                        flag until you have looked at your results without it, and
                        preferably tried to eliminate long genes by using better search
                        models.""")

    return parser

####### File conversion ########
def phylipise(species, number, size=8):
    padding = size - len(species) - len(str(number))
    if padding < 1:
        unpaddedname = "%s%s" % (species, number)
        shortname = unpaddedname[:10]
    else:
        shortname = "%s%s%s" % (species, "0" * padding, number)
    return shortname

def make_phylip(fastaalignment, logfile):
    "Convert a fasta file alignment to phylip format"
    phylip_alignment = logfile[:-3] + 'phylip'

    input_handle = open(fastaalignment, 'rb')
    output_handle = open(phylip_alignment, 'w')

    alignment = AlignIO.read( input_handle, "fasta")
    AlignIO.write(alignment, output_handle, "phylip")

    input_handle.close()
    output_handle.close()

    return phylip_alignment

def trim_name_dross(genename):
    if genename.find('|') >= 0:
        return genename[genename.find('|')+1:]
    else:
        return genename

####### fasta file operations ########
def get_gene_fastas(genes=None, species=None, fastafile=None,
                    specieslist = [], comment=None, short=False):
    """
    Can either be given as a transcript name to be searched within the peptide databases,
    or can be a fasta file.
    """
    if genes:
        for gene in genes:
            if species in specieslist:
                reportedspecies = species
                seq = gm.extractseq(gene, db=dbpaths[species + '_lpep'])

                if len(seq) == 0:
                    verbalise("R", "Transcript %s could not be extracted from the LNRP database for species %s" % (gene, species))
                    exit()
            else:   # if no species is given, check all LNRP files
                for sp in specieslist:
                    seq = gm.extractseq(gene, db=dbpaths[sp + '_lpep'])
                    if len(seq) > 0:   # found a match!
                        reportedspecies = sp
                        break
                else:
                    verbalise("R", "Transcript %s could not be extracted from the LNRP database for species %s" % (gene, species))
                    defline, seq, reportedspecies = None, None, None

            # create fasta file from extracted sequence:
            if short:
                name = phylipise(reportedspecies, short)
                defline = ">%s (%s) %s" % (name, reportedspecies, comment)
            else:
                defline = ">%s (%s) %s" % (gene, reportedspecies, comment)

            yield defline, seq, reportedspecies

    elif fastafile:
        handle = open(fastafile, 'rb')
        seq = ""
        for line in handle:
            if line[0] == '>':
                if seq != "":
                    yield defline, seq, None
                defline = line.strip()
                seq = ""
            else:
                seq += line.strip()
        else:
            yield defline, seq, None

    else:
        yield None, None, None

def get_similar_sequences(temp_dir, buildhmmer=False, fastafile=None,
                        specieslist={}, species=None,
                        mincollect=2, globalthresh=0.2, localthresh=0.8):
    if buildhmmer:
        hmminput = os.path.join(temp_dir, "hmminput.fa")
        handle = open(hmminput, 'w')
        seqcount = 0
        for defline, seq, species in get_gene_fastas(genes=genes,
                                                    species=None,
                                                    fastafile=fastafile,
                                                    specieslist=specieslist):
            seqcount += 1
            fasta_seq = "%s\n%s\n" % (defline, seq)
            handle.write(fasta_seq)
        handle.close()

        # create alignment of input sequences:
        verbalise("B",
                "Creating alignment and hmmer model of %d input sequences" % seqcount)
        mafft_align1 = os.path.join(temp_dir, "mafft_align_input.fa")
        mafft_align(hmminput, mafft_align1)

        # create hmmbuild model of alignment:
        hmmmodel = os.path.join(temp_dir, "hmmmodel.fa")
        open(hmmmodel, 'a').close()
        handle = os.popen(" ".join(['hmmbuild --informat afa', hmmmodel, mafft_align1]))
        handle.close()

        homologlist = hmmer_search(None,
                                    specieslist,
                                    query_species=species,
                                    minthresh=localthresh,
                                    temp_dir=temp_dir,
                                    mincollect=mincollect,
                                    globalthresh=globalthresh,
                                    hmmfile=hmmmodel)

        os.remove(mafft_align1)
        os.remove(hmminput)

    else:
        # run phmmer on a single input gene/sequence:
        for defline, seq, species in get_gene_fastas(genes=genes,
                                                    species=species,
                                                    fastafile=fastafile,
                                                    specieslist=specieslist):
            fasta_seq = "%s\n%s\n" % (defline, seq)
            verbalise("C", fasta_seq)

        ## phmmer all lpep files
        homologlist = hmmer_search(fasta_seq,
                                    specieslist,
                                    query_species=species,
                                    minthresh=localthresh,
                                    temp_dir=temp_dir,
                                    mincollect=mincollect,
                                    globalthresh=globalthresh,
                                    hmmfile=None)

    return homologlist

def mafft_align(inputfasta, outputfasta):
    mafft = os.popen( 'mafft --quiet ' + inputfasta )
    handle = open( outputfasta, 'w')
    for line in mafft:
        handle.write(line)
    handle.close()
    mafft.close()

def rank_scores(homologlist, thresh1, thresh2=None, genename=None, outfile=None, showplot=False):
    yvalues = sorted([val[1] for val in homologlist.values()], reverse=True)
    plt.plot(yvalues)
    score_cutoff = thresh1 * max(yvalues)
    sample_cutoff = sum(1 for s in yvalues if s >= thresh1 * max(yvalues))
    plt.axhline( score_cutoff , color='r' )
    if thresh2:
        plt.axhline( thresh2 * max(yvalues) , color='r' )
    plt.axvline( sample_cutoff -1 , color='g' )
    plt.text(sample_cutoff + 1,score_cutoff + 10 , "(%d,%d)" % (sample_cutoff,score_cutoff) )
    plt.xlabel("Gene rank")
    plt.ylabel("Phmmer score")
    plt.title("Ranking of phmmer scores for alignment with %s" % genename)
    if outfile:
        plt.savefig(outfile, format='png')
    if showplot:
        plt.show()
    else:
        plt.close()

def find_holes(seq):
    allholes = re.findall('([A-Za-z]+)(-{50,}[A-Za-z])', seq)
    dists = []
    for hole in allholes:
        dists.append(len(hole[0]))
        dists.append(len(hole[1]))
    return dists

def find_biggest_hole(seq):
    allholes = re.findall('-+', seq)
    if len(allholes) > 0:
        biggesthole = len(sorted(allholes)[-1])
        pattern = '(.+)-{' + str(biggesthole) + '}(.+)'
        bigsearch = re.search(pattern, seq)
        return len(bigsearch.group(1)), biggesthole, len(bigsearch.group(2))
    else:
        return 0,0,len(seq)

def get_pcmatch(seq):
    if len(seq) == 0:
        return 0, 0
    minigaps = len(re.findall('-', seq))
    width = len(seq)
    matches = width - minigaps
    assert matches >= 0
    pcmatch = 10.0 * matches / width
    return width, pcmatch

def display_alignment(fastafile, conversiondic={}, outfile=None, showplot=True):
    """
    Draw an alignment graph in the vein of BLAST alignment results on NCBI.
    colour scale represents the % match as base 10, to allow flooring of actual percentages
    to find appropriate colour group. Key number represents the minimum value the % match
    must be to join that colour group.
    """
    cmap = {0:'white', 1:'silver', 2:'tan' ,3:'cornflowerblue' ,4:'blue' ,5:'darkcyan' ,6:'green', 7:'gold' ,8:'orangered' ,9:'red' ,10:'maroon'}


    graph_points = {}
    hole_points = {}
    for defline, seq, species in get_gene_fastas(fastafile=fastafile):
        # determine the smallest reportable gap size is:
        repgap = int(len(seq)/10)

        # get distances and coverage percentages
        points = re.search('^(-*)(\S+[A-Za-z])(-*)$', seq)
        if points:
            pattern = '-{' + str(repgap) + ',}'
            fragments = re.findall(pattern, points.group(2))
            # set starting pos to beginning of matching sequence:
            spos = len(points.group(1))
            if len(fragments) > 0:
                """
                dists is a list of tuples, each tuple containing the start position of
                a large gap,  the length of the gap, the start of the preceding non-gap
                fragment, its width and the % match.
                """
                dists = []
                for frag in fragments:
                    nextgap = seq.find(frag, spos)
                    width, pcmatch = get_pcmatch(seq[spos:nextgap])
                    dists.append((nextgap,len(frag), spos, width, pcmatch))
                    spos = nextgap + len(frag)

                else:
                    lastfrag = points.group(3)
                    nextgap = len(seq) - len(lastfrag)
                    width, pcmatch = get_pcmatch(seq[spos:nextgap])
                    dists.append((0,0, spos, width, pcmatch))

            else:
                width, pcmatch = get_pcmatch(points.group(2))
                dists = [(0,0,spos, width, pcmatch)]

        else:
            dists = [(0,0,0,1,0)]

        # get name (convert if possible):
        namesearch = re.search('>(\S+)', defline)
        if namesearch:
            genename = namesearch.group(1)
        else:
            genename = defline.strip()
        if genename in conversiondic:
            fullname = conversiondic[genename][0]
        else:
            fullname = genename

        graph_points[fullname] = dists

    # get coords for alignment:
    keynames = sorted(graph_points.keys(), reverse=True)
    name_pos = np.arange(len(graph_points)) + 0.5
    y_frame = { k:y for k,y in zip(keynames, name_pos)}

    y_pos, lefts, widths, colors, bh_lefts, bh_widths = [], [], [], [], [], []
    for k in keynames:
        for dists in graph_points[k]:
            y_pos.append(y_frame[k])
            lefts.append(dists[2])
            widths.append(dists[3])
            colors.append(cmap[floor(dists[4])])
            bh_lefts.append(dists[0])
            bh_widths.append(dists[1])

    # plot graph:
    if 75 > len(keynames) > 25:
        plt.figure(figsize=(10,10))
    elif len(keynames) >= 75:
        plt.figure(figsize=(10,20))
    plt.barh(left=lefts,    width=widths,    bottom=y_pos, height=0.8, color=colors)
    plt.barh(left=bh_lefts, width=bh_widths, bottom=y_pos, height=0.8, color='white',
            alpha=0.5)
    plt.yticks(name_pos + 0.4, keynames)
    plt.xlabel("position (aa)")
    plt.title("Alignment of genes")
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, format='png')
    if showplot:
        plt.show()
    else:
        plt.close()

####### hmmer functions ########
def hmmer_search(fasta_seq, specieslist, query_species,  temp_dir,
                    minthresh=0.8, mincollect=2, globalthresh=0.01, hmmfile=None):
    """
    HMMER search of longest non-redundant peptide fasta files using either phmmer (with
    query protein as input) or hmmersearch (using constructed hmmfile as input).
    Finds best score, and collects all proteins with score > PC% (default 80%) of the
    best score, but limited to no more than X (default 2) proteins per species.
    """

    # set search type and input files:
    if hmmfile:
        searchcmd = 'hmmsearch'
        hmminput = hmmfile
    else:
        searchcmd = 'phmmer'
        hmminput = os.path.join(temp_dir,"seq.fasta")
        handle = open(hmminput, 'w')
        handle.write(fasta_seq)
        handle.close()

    verbalise("B", "Finding homologs using %s..." % searchcmd)
    homologlist = {}
    all_results = {}
    filtered_results = {}
    has_bestscore = False
    for sp in [query_species] + [ s for s in specieslist if s != query_species ]:
        try:
            phandle = os.popen( " ".join([searchcmd, hmminput, lpep_paths[sp]]) )
        except KeyError:
            continue

        # parse phmmer results:
        all_results[sp] = parse_the_hmmer(phandle)

        # determine cutoff threshold for filtering results:
        if sp == query_species:
            bestscore = max( v[1] for v in all_results[sp].values())
            has_bestscore = True
        elif has_bestscore:
            pass
        else:
            bestscore = max( v[1] for v in all_results[sp].values())
        cutoff_thresh = minthresh * bestscore

        # filter local file based on parameters given:
        ars = all_results[sp]
        filtered_results[sp] = { gene:(sp, ars[gene][1]) for i,gene in enumerate(ars) if ars[gene][1] >= cutoff_thresh or i < mincollect}

    # filter for global threshold (ie, based on % of best match). Most useful if no
    # species has been specified for fasta file.
    bestscore = max( v[1] for s in filtered_results for v in filtered_results[s].values())
    global_thresh = globalthresh * bestscore
    verbalise("C", "Best hmmer score = %d" % bestscore)
    #verbalise("M", filtered_results.items())
    homologlist = { k:v for nd in filtered_results.values() for (k,v) in nd.items() if v[1] >= global_thresh }

    # clean up temporary files
    os.remove(hmminput)

    return homologlist

def parse_the_hmmer(handle):
    """
    parses the protein matches from a hmmer search and returns a dictionary of peptides
    and their associated score and p-value.
    """
    parse_dic = {}
    lcount = 0
    collected = 0
    for line in handle:
        lcount += 1
        if len(line) < 2:
            continue
        if line.split()[1] in ['hits', 'inclusion', 'annotation']:
            break
        else:
            try:
                score = float(line.split()[1])
                pvalue = eval(line.split()[0])
            except ValueError:
                continue
            else:
                parse_dic[line.split()[8]] = (pvalue, score)

    handle.close()
    return parse_dic

####### Phylogeny creation/manipulation ########
def raxml_phylogeny(phylip_alignment, logfile, bootstrap=False):
    if bootstrap:
        bootstrapopt = '-N 100 -x ' + str(bootstrap)
        bs_str = 'bs.'
    else:
        bootstrapopt = '-N 1'
        bs_str = ""

    # determine path of final RAxML output file
    raxml_outfile = os.path.basename(logfile[:-3] + bs_str + 'raxml.out')
    if bootstrap:
        prefix = "RAxML_bootstrap."
    else:
        prefix = "RAxML_bestTree."
    raxml_final = os.path.join(os.path.dirname(logfile), prefix + raxml_outfile)

    cmd = " ".join(["raxmlHPC", "-s", phylip_alignment,
                                    '-w', os.path.dirname(logfile),
                                    "-n", raxml_outfile,
                                    '-p', "12345",
                                    "-m", 'PROTGAMMALG',
                                    bootstrapopt,
                                    ])
    os.system( cmd )
    return raxml_final

def rename_newick(raxml_final, conversiondic={}):
    #replace short names in newick file with full names
    if os.path.exists(raxml_final):
        handle = open(raxml_final, 'rb')
        newfile = raxml_final[:-3] + "scored.nwk"
        newhandle = open(newfile, 'w')
        for line in handle: # should only be one line in file
            for shortname in conversiondic:
                line = re.sub(shortname,
                            conversiondic[shortname][0] + "_" + str(conversiondic[shortname][1]),
                            line )
            newhandle.write(line)
        handle.close()
        newhandle.close()
    else:
        newfile = None
    return newfile

def apply_boostrap(besttree, bootstraps, logfile):
    """
    Takes the bootstrap trees and the best tree and puts the bootstrap support onto the
    best tree.
    """
    outfile = os.path.basename(logfile[:-3] + 'raxml_merged.nwk')
    cmd = " ".join(["raxmlHPC",     '-z', bootstraps,
                                    '-w', os.path.dirname(logfile),
                                    '-n', outfile,
                                    '-m', 'PROTGAMMALG',
                                    '-t', besttree,
                                    '-f', 'b',
                                    ])
    os.system(cmd)
    return os.path.join(os.path.dirname(logfile), "RAxML_bipartitions." + outfile)

################################################

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    temp_dir = tempfile.mkdtemp()

    ############ if generating phylogeny only ###########
    if args.raxml_only:
        if not args.fasta:
            verbalise("R", "No fasta alignment was supplied!")
            exit()
        verbalise("B", "Running RAxML analysis to construct phylogeny using supplied alignment")
        phylip_alignment = make_phylip(args.fasta, logfile)
        raxml_final = raxml_phylogeny(phylip_alignment, logfile, bootstrap=args.bootstrap)
        if args.name_conversion:
            handle = open(args.name_conversion, 'rb')
            conv_dic = { line.split()[0]:(line.split()[2], line.split()[1]) for line in handle }
            handle.close()
            rename_newick(raxml_final, conversiondic=conv_dic)
        exit()

    # initialise dictionary of all accessible longest non-redundant peptide fasta files
    specieslist = [ "Ador", "Aech", "Aflo", "Amel", "Apis", "Aros", "Bimp", "Bmor",
                "Bter", "Cele", "Cflo", "Csol", "Dcit", "Fari", "Hsal", "Lhum",
                "Mdem", "Mpha", "Mrot", "Nvit", "Oabi", "Pbar", "Pcan", "Sinv",
                "Tcas", "Waur", "Cbir", "Ebur", "Dmel" ]
    lpep_paths = { s:dbpaths[s+'_lpep'] for s in specieslist }

    ######### Get protein sequences #########
    genes = config.make_a_list(args.gene)
    homologlist = get_similar_sequences(temp_dir,
                                        buildhmmer=args.buildhmmer,
                                        fastafile=args.fasta,
                                        specieslist=specieslist,
                                        species=args.species,
                                        mincollect=args.mincollect,
                                        globalthresh=args.globalthresh,
                                        localthresh=args.scorethresh)

    ######### Extract identified sequences from LNRP fasta files #########
    conv_handle = open(logfile[:-3] + 'name_conversion.txt', 'w')
    conv_dic = {}
    itercount = 0
    previousseq = ""
    seqdic = {}         # loaded up to remove duplicate sequences
    excluded_genes = config.make_a_list(args.exclude_genes)
    excluded_species = config.make_a_list(args.exclude_species)

    for homolog in sorted(homologlist):
        # remove excluded genes before bothering to look up their sequence:
        searchname = trim_name_dross(homolog)
        if searchname in excluded_genes:
            continue
        if homologlist[homolog][0] in excluded_species:
            continue

        # extract sequences of remaining genes and add to conversion dictionary
        itercount += 1

        for defline, seq, spec in get_gene_fastas(genes=[searchname],
                                    species=homologlist[homolog][0],
                                    fastafile=None,
                                    specieslist = specieslist,
                                    comment=homologlist[homolog][1],
                                    short=itercount):

            if args.maxlength and len(seq) > args.maxlength:
                continue
            seqdic[seq] = defline

        shortname = phylipise(homologlist[homolog][0], itercount)
        conv_handle.write("%s %-5d %s\n" % (shortname,
                                          homologlist[homolog][1],
                                          homolog))
        conv_dic[shortname] = (homolog, homologlist[homolog][1])
    conv_handle.close()

    verbalise("G", "%d non-excluded homologous sequences found" % len(seqdic))

    # write multifasta file:
    homolog_fasta = os.path.join(temp_dir,"homolog.fasta")
    handle = open(homolog_fasta, 'w')
    for seq in seqdic:
        fastaseq = "%s\n%s\n" % (seqdic[seq], seq)
        handle.write(fastaseq)
    handle.close()

    # show score curve for identified genes:
    rank_scores(homologlist, args.scorethresh, args.globalthresh, args.gene,
                logfile[:-3] + "ranking.png", showplot=args.display_on)
    # show raw protein sizes:
    display_alignment(  homolog_fasta,
                        conversiondic=conv_dic,
                        outfile=logfile[:-3] + 'homologs.png',
                        showplot=args.display_on)

    ######### MAFFT alignment of extracted sequences #########
    mafft_alignment = logfile[:-3] + 'mafft.fa'
    mafft_align(homolog_fasta, mafft_alignment)
    display_alignment(mafft_alignment,
                        conversiondic=conv_dic,
                        outfile=logfile[:-3] + 'mafft.png',
                        showplot=args.display_on)

    ######### RaXML phylogenetic analysis of alignment #########
    """
    Using PROT-LG-GAMMA model, a single tree and no bootstrapping (though all of these
    could be setup to allow overriding for fringe case analyses).
    """
    verbalise("B", "Running RAxML analysis to construct phylogeny")
    phylip_alignment = make_phylip(mafft_alignment, logfile)
    if args.bootstrap:
                raxml_best = raxml_phylogeny(phylip_alignment,
                                        logfile,
                                        bootstrap=False)
                best_renamed = rename_newick(raxml_best, conversiondic=conv_dic)
                raxml_bstrap = raxml_phylogeny(phylip_alignment,
                                        logfile,
                                        bootstrap=args.bootstrap)
                bstrap_renamed = rename_newick(raxml_bstrap, conversiondic=conv_dic)
                final_tree = apply_boostrap(best_renamed, bstrap_renamed, logfile)
                verbalise("Y",
                    "Best tree with bootstrap support can be found at %s" % final_tree)
    else:
        raxml_final = raxml_phylogeny(phylip_alignment, logfile, bootstrap=args.bootstrap)
        raxml_renamed = rename_newick(raxml_final, conversiondic=conv_dic)

    # clean up temp files and directory
    for file in [ homolog_fasta, ]:
        if os.path.exists(file):
            os.remove(file)

    os.rmdir(temp_dir)  # dir must be empty!





