#!/usr/bin/env python

"""
This is a script to calculate the length of each isoform of a gene from the gff file,
extract the peptide sequences from a pep file, then create a new pep file that only
contains the longest isoform from each gene.
"""

import sys
import re
import argparse

from genomepy import config

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "Performs set analysis on multiple lists and outputs venn diagrams.")

    # input options
    parser.add_argument("-g", "--gff", type=str,
                        help="specify the gff file to use")
    parser.add_argument("-f", "--fasta", type=str,
                        help="specify the fasta file to use")
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="don't print messages")
    parser.add_argument("-k", "--speckey", type=str, default="",
                        help="specify the species prefix")
    parser.add_argument("-o", "--output", type=str, default='lorf2.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-L", "--find_longest", action='store_true',
                        help="find longest peptide isoform for each gene in specified gff")
    parser.add_argument("-C", "--convert", action='store_true',
                        help="convert specified gff to gtf format")
    parser.add_argument("-A", "--audit", action='store_true',
                        help="count number of features in gff file")
    parser.add_argument("-R", "--replace", action='store_true',
                        help="replace isoforms in fasta def line with genes from gff file")
    parser.add_argument("-M", "--replace_all", action='store_true',
                        help="replace isoforms in 'fasta' columns with genes from gff file")
    parser.add_argument("-b", "--large_term", type=str, default='gene',
                        help = "the term to replace to")
    parser.add_argument("-s", "--small_term", type=str, default='protein_id',
                        help = "the term to replace from")



    return parser

def build_isodic(gff, idbig='gene', idsmall='protein_id'):
    handle = open(gff, 'rb')
    genedic = {}
    for line in handle:
        attributes = attr(line)
        genedic = update_genedic(genedic, attributes, True, idbig=idbig, idsmall=idsmall)
    handle.close()
    return genedic

def replace_all(isodic, fasta, outfile, speckey):
    handle = open(fasta, 'rb')
    outhandle = open(outfile, 'w')

    for line in handle:
        # check the line has the species prefix:
        columns = line.split()
        specheck = [ re.search(speckey, col) for col in columns]
        pattern = speckey + "\|?([\S]+)"
        isocheck = None
        if specheck[0]:
            isocheck = re.search(pattern, columns[0])
            col = 0
        elif  specheck[1]:
            isocheck = re.search(pattern, columns[1])
            col = 1

        if isocheck:
            if isocheck.group(1) in isodic:
                columns[col] = speckey + "|" + isodic[isocheck.group(1)][0]
                try:
                    outhandle.write("%s\n" % ("\t".join(columns)))
                except TypeError:
                    verbalise("R", columns)
                    exit()
            else:
                outhandle.write(line)
        else:
            outhandle.write(line)

def replace_fasta(isodic, fasta, outfile):
    handle = open(fasta, 'rb')
    outhandle = open(outfile, 'w')

    for line in handle:
        if line[0] == '>':
            predef = line.split()[0][:6]
            isodef = line.split()[0][6:]
            if isodef in isodic:
                outhandle.write("%s%s\n" % (predef, isodic[isodef][0]))
            else:
                outhandle.write(line)
        else:
            outhandle.write(line)
    outhandle.close()
    handle.close()

def findlongest(gff, pepfile):
    isolen = {}

    gff_h = open(gff, 'rb')
    for line in gff_h:
        if line[0]=='#':
            continue
        if line.split()[2] == 'CDS':
            cdstart = int(line.split()[3])
            cdstop = int(line.split()[4])
            attr = { term.split('=')[0]:term.split('=')[1] for term in " ".join(line.split()[8:]).split(';')  }
            try:
                gene = attr['gene']
                mrna = attr['protein_id']
            except KeyError:
                if re.search('RNA|CDS|gene',line.split()[2]):
                    print line
                    print attr
                    print "-----"
                continue
            if gene not in isolen:
                isolen[gene] = { mrna:abs(cdstart-cdstop) }
            elif mrna not in isolen[gene]:
                isolen[gene][mrna] = abs(cdstart-cdstop)
            else:
                isolen[gene][mrna] += abs(cdstart-cdstop)

    gff_h.close()


    # load peptide sequences:
    print "Loading all peptide sequences."

    pep_h = open(pepfile, 'rb')
    pepseq = {}
    peptide = 'Should not be here'
    for line in pep_h:
        if line[0] == '>':
            try:
                peptide = re.search('(gb|ref)\|([^\|]+)\|', line).group(2)
            except AttributeError:
                print line
                peptide = "error"
            if peptide == None:
                print line
                print "Peptide read as NoneType!!!"
            pepseq[peptide] = ""
        else:
            pepseq[peptide] += line.strip()



    # write longest peptide sequences to fasta file:
    longest = gff[:-3]+"longest.pep"
    print "writing longests sequences to ", longest
    longest_h = open(longest, 'w')
    for gene in isolen:
        for isoform in isolen[gene]:
            if isolen[gene][isoform] == max(isolen[gene].values()):
                longest_h.write(">Cbir|%s\n" % (isoform))
                try:
                    longest_h.write(pepseq[isoform] + "\n")
                except KeyError:
                    print isoform
                    continue
                break

    longest_h.close()

def attr(line):
    if line[0] == '#':
        return None
    cols = line.split()
    attr = { pair.split('=')[0]:pair.split('=')[1] for pair in " ".join(cols[8:]).split(';')  }
    return attr

def update_genedic(genedic, attr, reverse=False, idbig='gene', idsmall='protein_id'):
    try:
        prot = attr[idsmall]
        gene = attr[idbig]
    except:
        return genedic
    else:
        if reverse:
            prot = attr[idbig]
            gene = attr[idsmall]
        if gene in genedic:
            if prot in genedic[gene]:
                return genedic
            else:
                genedic[gene] += [prot]
        else:
            genedic[gene] = [prot]
    return genedic

def convert_gff(gff, outfile):
    handle = open(gff, 'rb')
    outhandle = open(outfile, 'w')

    exoncount = 0
    writecount = 0
    for line in handle:
        if line[0] == '#':
            continue
        attributes = attr(line)
        if line.split()[2]  == 'exon':
            exoncount += 1
            try:
                outhandle.write('%s\tgene_id "%s";\ttranscript_id "%s"; id "%s";\n' % (
                                "\t".join(line.split()[:8]),
                                attributes['gene'],
                                attributes['transcript_id'],
                                attributes['ID'])
                            )
            except KeyError as ke:
                if attributes['gbkey'] == 'exon':
                    pass
                elif attributes['gbkey'] == 'tRNA':
                    outhandle.write('%s\tgene_id "%s";\ttranscript_id "%s"; id "%s";\n' % (
                                "\t".join(line.split()[:8]),
                                attributes['gene'],
                                attributes['Dbxref'],
                                attributes['ID'])
                            )
                    writecount += 1
                else:
                    verbalise("R", "KeyError for %s:\n" % (ke), attributes)
            else:
                writecount += 1
    outhandle.close()
    handle.close()
    return exoncount, writecount

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    if args.convert:
        exons, wexons = convert_gff(args.gff, logfile[:-3]+'gtf')
        verbalise("G", "%d/%d exons successfully written to file" % (wexons, exons))
    if args.audit:
        genedic = {}
        handle = open(args.gff, 'rb')
        for line in handle:
            newattr = attr(line)
            genedic = update_genedic(genedic, newattr, idbig=args.large_term, idsmall=args.small_term)

        verbalise( len(genedic), genedic.keys()[:5], genedic.values()[:5] )
        verbalise( sum( 1 for val in genedic if len(genedic[val]) > 1) )
    elif args.find_longest:
        findlongest(args.gff, logfile[:-3]+'out')
    elif args.replace:
        verbalise("building isoform dictionary")
        isodic = build_isodic(args.gff, idbig=args.large_term, idsmall=args.small_term)
        verbalise("G", "%d isoforms added" % (len(isodic)))
        egkeys = " : ".join(isodic.keys()[:5])
        egvals = " : ".join([str(val) for val in isodic.values()[:5]])
        verbalise("Y", "%s\n%s" % (egkeys, egvals) )
        verbalise("writing new file")
        replace_fasta(isodic, args.fasta, logfile[:-3]+'out')
    elif args.replace_all:
        verbalise("building isoform dictionary")
        isodic = build_isodic(args.gff, idbig=args.large_term, idsmall=args.small_term)
        verbalise("G", "%d isoforms added" % (len(isodic)))
        egkeys = " : ".join(isodic.keys()[:5])
        egvals = " : ".join([str(val) for val in isodic.values()[:5]])
        verbalise("Y", "%s\n%s" % (egkeys, egvals) )
        verbalise("writing new file")
        replace_all(isodic, args.fasta, logfile[:-3]+'out', args.speckey)
