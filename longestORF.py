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
            "Performs file modifications on gff and fasta files.")
    # logging options
    parser.add_argument("-o", "--output", type=str, default='lorf2.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")

    # input options
    parser.add_argument("-g", "--gff", type=str,
                        help="specify the gff file to use")
    parser.add_argument("-f", "--fasta", type=str,
                        help="specify the fasta file to use")

    # analysis options
    parser.add_argument("-p", "--pattern", type=str, default='(gb|ref)\|([^\|]+)\|',
                        help="specify the regex to locate peptide name")
    parser.add_argument("-c", "--column", type=int, default=0,
                        help = "the column in which to convert names. (default=0)")
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="don't print messages")
    parser.add_argument("-k", "--speckey", type=str, default="",
                        help="specify the species prefix")
    parser.add_argument("-L", "--find_longest", action='store_true',
                        help="find longest peptide isoform for each gene in specified gff")
    parser.add_argument("-C", "--convert", action='store_true',
                        help="convert specified gff to gtf format")
    parser.add_argument("-A", "--audit", action='store_true',
                        help="count number of features in gff file")
    parser.add_argument("-R", "--replace", action='store_true',
                        help="replace isoforms in fasta def line with genes from gff file")
    parser.add_argument("-M", "--replace_all", action='store_true',
                        help="replace isoforms in specified column (use --column and --fasta) with genes from gff file")
    parser.add_argument("-b", "--large_term", type=str, default='gene',
                        help = "the term to replace to. (Default='gene')")
    parser.add_argument("-s", "--small_term", type=str, default='protein_id',
                        help = "the term to replace from. (Default='protein_id'")
    parser.add_argument("-a", "--addname", type=str, default='gene',
                        help = """In the specified gff file, replace the large term in
                        the specified level with the small term""")



    return parser

def build_isodic(gff, idbig='gene', idsmall='protein_id', reverse=True):
    handle = open(gff, 'rb')
    genedic = {}
    for line in handle:
        attributes = attr(line)
        genedic = update_genedic(genedic, attributes, reverse=reverse, idbig=idbig, idsmall=idsmall)
    handle.close()
    return genedic

def col_counter(columns,column):
    if column:
        yield columns[column]
    else:
        for col in columns:
            yield col

def replace_all(isodic, filename, outfile, prefix, column=None):
    """
    search for gene ids in the specified columns (if none, will search all columns)



    """
    handle = open(filename, 'rb')
    outhandle = open(outfile, 'w')

    for line in handle:
        # check the line has the species prefix:
        columns = line.split()
        prefix_check = [ re.search(prefix, col) for col in columns]
        pattern = prefix + "\|?([\S]+)"
        isocheck = None

        valid_columns = col_counter(columns, column)
        for col, contents in enumerate(valid_columns):
            if prefix_check[col]:
                isocheck = re.search(pattern, columns[col])

                if isocheck:
                    if isocheck.group(1) in isodic:
                        columns[col] = prefix + "|" + isodic[isocheck.group(1)][0]
        try:
            outhandle.write("%s\n" % ("\t".join(columns)))
        except TypeError:
            verbalise("R", columns)
            exit()
        #            else:
        #                outhandle.write(line)
        #        else:
        #            outhandle.write(line)

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

def findlongest(gff, pepfile, species_key, pattern='(gb|ref)\|([^\|]+)\|'):
    print "Searching using pattern %r" % pattern
    isolen = {}
    err_cnt = 0

    gff_h = open(gff, 'rb')
    for line in gff_h:
        if line[0]=='#':
            continue
        if line.split()[2] == 'CDS':
            cdstart = int(line.split()[3])
            cdstop = int(line.split()[4])
            attr = { term.split('=')[0]:term.split('=')[1] for term in " ".join(line.split()[8:]).split(';')  }
            
            if 'gene' in attr:
                gene = attr['gene']
            elif 'Name' in attr:
                gene = attr['Name']
            else:
                err_cnt += 1
                continue
            
            if 'protein_id' in attr:
                mrna = attr['protein_id']
            elif 'transcript_id' in attr:
                mrna = attr['transcript_id']
            else:
                err_cnt += 1
                continue
            
            if gene not in isolen:
                isolen[gene] = { mrna:abs(cdstart-cdstop) }
            elif mrna not in isolen[gene]:
                isolen[gene][mrna] = abs(cdstart-cdstop)
            else:
                isolen[gene][mrna] += abs(cdstart-cdstop)

    if err_cnt > 0:
        verbalise("R", err_cnt, "CDS lines could not be parsed.")
    gff_h.close()
    verbalise("G", "%d genes parsed successfully" % len(isolen))
    print isolen.keys()[0], isolen[isolen.keys()[0]]

    # load peptide sequences:
    verbalise("B", "Loading all peptide sequences.")

    pep_h = open(pepfile, 'rb')
    pepseq = {}
    peptide = 'Should not be here'
    for line in pep_h:
        if line[0] == '>':
            # identify peptide name: default is NCBI format >gi|328776094|ref|XP_001122629.2| PRED..
            try:
                peptide = re.search(pattern, line).group(2)
            except AttributeError:
                print line.strip()
                peptide = "error"
            if peptide == None:
                print line
                print "Peptide read as NoneType!!!"
            pepseq[peptide] = ""
        else:
            pepseq[peptide] += line.strip()
    pep_h.close()
    verbalise("G", "%d peptide sequences successfully extracted" % len(pepseq))
    print pepseq.keys()[:5]

    # write longest peptide sequences to fasta file:
    longest = gff[:-4]+"longest.pep"
    verbalise("Y", "writing longests sequences to ", longest)
    longest_h = open(longest, 'w')
    for gene in isolen:
        for isoform in isolen[gene]:
            if isolen[gene][isoform] == max(isolen[gene].values()):
                longest_h.write(">%s|%s\n" % (species_key, isoform))
                try:
                    longest_h.write(pepseq[isoform] + "\n")
                except KeyError:
                    print isoform
                    continue
                break

    longest_h.close()

def attr(line):
    if len(line) == 0:
        return None
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
                    try:
                        outhandle.write('%s\tgene_id "%s";\ttranscript_id "%s"; id "%s";\n' % (
                                    "\t".join(line.split()[:8]),
                                    attributes['gene'],
                                    attributes['Dbxref'],
                                    attributes['ID'])
                                )
                    except:
                        print line
                    else:
                        writecount += 1
                else:
                    verbalise("R", "KeyError for %s:\n" % (ke), attributes)
            else:
                writecount += 1
    outhandle.close()
    handle.close()
    return exoncount, writecount

def replace_gff(isodic, gff_file, outfile,
                    level='gene', idbig='gene', idsmall='product'):

    "replaces idbig name in level with idsmall"
    handle = open(gff_file, 'rb')
    outhandle = open(outfile, 'w')

    foundidbig = []

    for line in handle:
        cols = line.split()
        if cols[2] == level:
            atts = attr(line)
            if idbig in atts:
                if atts[idbig] in isodic:
                    atts[idbig] = isodic[atts[idbig]][0]
                    try:
                        cols[8] = ";".join([ k + "=" + v for k,v in atts.items() ])
                    except TypeError:
                        verbalise("R", atts.items())
                elif 'pseudo' in atts:
                    atts[idbig] = "pseudo_" + atts[idbig]
                    try:
                        cols[8] = ";".join([ k + "=" + v for k,v in atts.items() ])
                    except TypeError:
                        verbalise("R", atts.items())
                else:
                    foundidbig.append(atts[idbig])
            else:
                cols[8] = " ".join(cols[8:])
        else:
            cols[8] = " ".join(cols[8:])
        outhandle.write("\t".join(cols[:9]) + "\n")
    handle.close()
    outhandle.close()

    if len(foundidbig) > 0:
        verbalise("R", "%d genes were not converted! (%s...)" % (len(foundidbig),",".join(foundidbig[:3])))


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
        findlongest(args.gff, args.fasta, args.speckey, pattern=args.pattern)
    elif args.replace:
        verbalise("building isoform dictionary")
        isodic = build_isodic(args.gff, idbig=args.large_term, idsmall=args.small_term)
        verbalise("G", "%d isoforms added" % (len(isodic)))
        egkeys = " : ".join(isodic.keys()[:5])
        egvals = " : ".join([str(val) for val in isodic.values()[:5]])
        verbalise("B", "%s\n%s" % (egkeys, egvals) )
        verbalise("writing new file")
        replace_fasta(isodic, args.fasta, logfile[:-3]+'out')
    elif args.replace_all:
        verbalise("building isoform dictionary")
        isodic = build_isodic(args.gff, idbig=args.large_term, idsmall=args.small_term)
        verbalise("G", "%d isoforms added" % (len(isodic)))
        egkeys = " : ".join(isodic.keys()[:5])
        egvals = " : ".join([str(val) for val in isodic.values()[:5]])
        verbalise("B", "%s\n%s" % (egkeys, egvals) )
        verbalise("writing new file")
        replace_all(isodic, args.fasta, logfile[:-3]+'out', args.speckey)
    elif args.addname:
        verbalise("building isoform dictionary")
        isodic = build_isodic(args.gff, idbig=args.large_term, idsmall=args.small_term,
                                reverse=False)
        verbalise("G", "%d isoforms added" % (len(isodic)))
        egkeys = " : ".join(isodic.keys()[:5])
        egvals = " : ".join([str(val) for val in isodic.values()[:5]])
        verbalise("B", "%s\n%s" % (egkeys, egvals) )
        replace_gff(isodic, args.gff, logfile[:-3]+'out',
                    idbig=args.large_term, idsmall=args.small_term,
                    level=args.addname)