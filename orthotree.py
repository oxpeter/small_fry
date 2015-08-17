#! /usr/bin/env python
"""
A wrapper to identify, extract, align and phylogenise a target gene. The primary purpose
of this module is for a user to confirm orthology of their gene, not to understand in
detail the phylogenetic relationships between all closely related proteins.
"""

import os
import sys
import tempfile
import argparse
import linecache

from Bio import AlignIO
import matplotlib.pyplot as plt

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

    # analysis options:
    parser.add_argument("-c", "--mincollect", type=int, default=2,
                        help="""minimum number of homologs to collect from each species,
                        regardless of how poor the similarity""")
    parser.add_argument("-t", "--scorethresh", type=float, default=0.8,
                        help="""minimum fraction of the score of the best matching
                        gene for a species for including additional matches""")
    parser.add_argument("-p", "--threads", type=int, default=2,
                        help="number of threads to use for RAxML calculations")
    parser.add_argument("-b", "--bootstrap", type=int,
                        help="""Perform rapid bootstrap analysis. Requires a random
                         integer for the bootstrap seed.""")
    parser.add_argument("-a", "--globalthresh", type=int,
                        help="""Identify the top proportion of all results globally.""")


    return parser

def phylipise(species, number, size=8):
    padding = size - len(species) - len(str(number))
    if padding < 1:
        unpaddedname = "%s%s" % (species, number)
        shortname = unpaddedname[:10]
    else:
        shortname = "%s%s%s" % (species, "0" * padding, number)
    return shortname

def trim_name_dross(genename):
    if genename.find('|') >= 0:
        return genename[genename.find('|')+1:]
    else:
        return genename

def get_gene_fasta(gene=None, species=None, fastafile=None,
                    specieslist = [], comment=None, short=False):
    """
    Can either be given as a transcript name to be searched within the peptide databases,
    or can be a fasta file.
    """
    if gene:
        if species in specieslist:
            seq = gm.extractseq(gene, db=dbpaths[species + '_lpep'])

            if len(seq) == 0:
                verbalise("R", "Transcript %s could not be extracted from the LNRP database for species %s" % (gene, species))
                exit()
        else:   # if no species is given, check all LNRP files
            for sp in specieslist:
                seq = gm.extractseq(gene, db=dbpaths[sp + '_lpep'])
                if len(seq) > 0:   # found a match!
                    species = sp
                    break
            else:
                verbalise("R", "Transcript %s could not be extracted from the LNRP database for species %s" % (gene, species))
                exit()

        # create fasta file from extracted sequence:
        if short:
            name = phylipise(species, short)
            defline = ">%s (%s) %s" % (name, species, comment)
        else:
            defline = ">%s (%s) %s" % (gene, species, comment)

    elif fastafile:
        handle = open(fastafile, 'rb')
        count = 0
        seq = ""
        for line in handle:
            if count == 0:
                defline = line.strip()
            else:
                seq += line.strip()

    else:
        defline = None
        seq = None

    return defline, seq, species

def phmmer_search(fasta_seq, specieslist, query_species,  temp_dir, minthresh=0.8, mincollect=2):
    """
    PHMMER of longest non-redundant peptide fasta files using query protein
    Find best score, and collect all proteins at with score > PC% (default 80%) of the
    best score, but limited to no more than X (default 5) proteins per species.
    """


    fasta_input = os.path.join(temp_dir,"seq.fasta")
    handle = open(fasta_input, 'w')
    handle.write(fasta_seq)
    handle.close()

    verbalise("B", "Finding homologs using phmmer...")
    homologlist = {}
    bestscore = None
    for sp in [query_species] + [ s for s in specieslist if s != query_species ]:
            phandle = os.popen( " ".join(['phmmer', fasta_input, lpep_paths[sp]]) )

            # parse output
            lcount = 0
            collected = 0
            for line in phandle:
                lcount += 1
                if len(line) < 2:
                    continue
                if sp == query_species and not bestscore:
                    try:
                        bestscore = float(line.split()[1])
                    except ValueError:
                        continue
                    verbalise("C", "Query alignment to own peptide:\n%s" % line.strip())
                    homologlist[line.split()[8]] = (sp, bestscore)
                else:
                    if line.split()[1] in ['hits', 'inclusion', 'annotation']:
                        break
                    else:
                        try:
                            currentscore = float(line.split()[1])
                        except ValueError:
                            continue

                    if currentscore >= minthresh * bestscore or collected < mincollect:
                        homologlist[line.split()[8]] = (sp, currentscore)
                        collected += 1
                    else:
                        break
            phandle.close()


    # clean up temporary files
    os.remove(fasta_input)

    return homologlist

def rank_scores(homologlist, thresh, genename):
    yvalues = sorted([val[1] for val in homologlist.values()], reverse=True)
    plt.plot(yvalues)
    score_cutoff = thresh * max(yvalues)
    sample_cutoff = sum(1 for s in yvalues if s >= thresh * max(yvalues))
    plt.axhline( score_cutoff , color='r' )
    plt.axvline( sample_cutoff -1 , color='g' )
    plt.text(sample_cutoff + 1,score_cutoff + 10 , "(%d,%d)" % (sample_cutoff,score_cutoff) )
    plt.xlabel("Gene rank")
    plt.ylabel("Phmmer score")
    plt.title("Ranking of phmmer scores for alignment with %s" % genename)
    plt.show()

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    temp_dir = tempfile.mkdtemp()

    # initialise dictionary of all accessible longest non-redundant peptide fasta files
    specieslist = [ "Ador", "Aech", "Aflo", "Amel", "Apis", "Aros", "Bimp", "Bmor",
                "Bter", "Cele", "Cflo", "Csol", "Dcit", "Fari", "Hsal", "Lhum",
                "Mdem", "Mpha", "Mrot", "Nvit", "Oabi", "Pbar", "Pcan", "Sinv",
                "Tcas", "Waur", "Cbir", "Ebur", "Dmel" ]
    lpep_paths = { s:dbpaths[s+'_lpep'] for s in specieslist }


    ######### Get protein ID #########
    # TODO: add ability to find transcript from gene name, to prevent having to guess the
    # appropriate transcript ID.
    defline, seq, species = get_gene_fasta(gene=args.gene,
                                        species=args.species,
                                        fastafile=args.fasta,
                                        specieslist=specieslist)
    fasta_seq = "%s\n%s\n" % (defline, seq)
    verbalise("C", fasta_seq)


    ## phmmer all lpep files
    homologlist = phmmer_search(fasta_seq,
                                specieslist,
                                query_species=species,
                                minthresh=args.scorethresh,
                                temp_dir=temp_dir,
                                mincollect=args.mincollect)



    ######### Extract identified sequences from LNRP fasta files #########
    conv_handle = open(logfile[:-3] + 'name_conversion.txt', 'w')
    itercount = 0
    previousseq = ""
    seqdic = {}  # loaded up to removed duplicate sequences

    if args.display_on:
        rank_scores(homologlist, args.scorethresh, args.gene)

    for homolog in sorted(homologlist):
        itercount += 1
        searchname = trim_name_dross(homolog)
        defline, seq, spec = get_gene_fasta(gene=searchname,
                                    species=homologlist[homolog][0],
                                    fastafile=None,
                                    specieslist = specieslist,
                                    comment=homologlist[homolog][1],
                                    short=itercount)
        seqdic[seq] = defline
        conv_handle.write("%s %-5d %s\n" % (phylipise(homologlist[homolog][0], itercount),
                                          homologlist[homolog][1],
                                          homolog))
    conv_handle.close()

    verbalise("G", "%d homologous sequences found" % len(seqdic))
    # write multifasta file:
    homolog_fasta = os.path.join(temp_dir,"homolog.fasta")
    handle = open(homolog_fasta, 'w')
    for seq in seqdic:
        fastaseq = "%s\n%s\n" % (seqdic[seq], seq)
        handle.write(fastaseq)
    handle.close()


    ######### MAFFT alignment of extracted sequences #########
    mafft = os.popen( 'mafft --quiet ' + homolog_fasta )
    mafft_alignment = logfile[:-3] + 'mafft.fa'
    phylip_alignment = logfile[:-3] + 'phylip'
    handle = open( mafft_alignment, 'w')
    for line in mafft:
        handle.write(line)
    handle.close()
    mafft.close()

    ######### RaXML phylogenetic analysis of alignment #########
    """
    Using PROT-LG-GAMMA model, a single tree and no bootstrapping (though all of these
    could be setup to allow overriding for fringe case analyses).
    """

    # convert fasta alignment to phylip
    input_handle = open(mafft_alignment, 'rb')
    output_handle = open(phylip_alignment, 'w')

    alignment = AlignIO.read( input_handle, "fasta")
    AlignIO.write(alignment, output_handle, "phylip")

    input_handle.close()
    output_handle.close()

    # run RAxML
    verbalise("B", "Running RAxML analysis to construct phylogeny")
    if args.bootstrap:
        bootstrapopt = '-x ' + str(args.bootstrap)
    else:
        bootstrapopt = ''

    cmd = " ".join(["raxmlHPC", "-s", phylip_alignment,
                                    '-w', os.path.dirname(logfile),
                                    "-n", os.path.basename(logfile[:-3] + 'raxml.out'),
                                    '-p', '54321',
                                    "-m", 'PROTGAMMALG',
                                    #'-T', str(args.threads),
                                    '-N', '1',
                                    bootstrapopt,
                                    ])
    verbalise("M", cmd)
    os.system( cmd )

    # clean up
    os.remove(homolog_fasta)
    os.rmdir(temp_dir)  # dir must be empty!





