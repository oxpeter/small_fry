#!/usr/bin/env python

"""
script to find the longest ORF from each transdecoder result from each EST from the
honey bee brain EST library BB16
"""

import re
import os, sys



def parse_name(line):
    # get details of next peptide:
    tag_s = re.search('cds\.([\w\.]*)\|', line)
    length_s = re.search('len:([0-9]*)', line)
    if tag_s is not None and length_s is not None:
        tag = tag_s.group(1)
        length = length_s.group(1)
    else:
        tag = "null_tag"
        length = 0
    return tag, length

def update_orfs(orfs, orf_sizes, tag, length, orf_string):
    # save details of previous peptide
    try:
        if length > orf_sizes[tag]:
            orf_sizes[tag] =  length
            orfs[tag] = orf_string
    except KeyError:
        orf_sizes[tag] = length
        orfs[tag] = orf_string
    return orfs, orf_sizes

def longest_pep(pepfile):
    orfs = {}
    orf_sizes = {}
    tag = 'null_tag'
    length = 0
    orf_string = ""

    pep_h = open(pepfile, 'rb')
    for line in pep_h:

        if line[0] == '>':
            # save details of previous peptide:
            orfs, orf_sizes = update_orfs(orfs, orf_sizes, tag, length, orf_string)
            # get values of new peptide:
            tag, length = parse_name(line)
            # reset some values:
            orf_string = ""

        else:
            # keep growing the peptide string:
            orf_string = orf_string + line
    else:
        orfs, orf_sizes = update_orfs(orfs, orf_sizes, tag, length, orf_string)

    return orfs

def write_to_file(orfs, outfile="apis.tag.fasta.transdecoder.longest.pep" ):

    # at end, write to file:
    outfile_h = open( outfile , 'w'
                )
    for tag in orfs:
        outfile_h.write(">%s\n%s\n" % (tag, orfs[tag]) )

def clean_defs(pepfile, outfile="apis.tag.fasta.transdecoder.clean.pep"):
    orfs = {}
    tag = 'null_tag'
    length = 0

    pep_h = open(pepfile, 'rb')
    outfile_h = open(outfile, 'w')
    for line in pep_h:
        if line[0] == '>':
            tag, length = parse_name(line)
            outfile_h.write(">%s\n" % (tag))
        else:
            outfile_h.write(line)

if __name__ == '__main__':
    if len(sys.argv) ==  3:
        pepfile = sys.argv[2]
    else:
        pepfile = 'apis.tag.fasta.transdecoder.pep'

    if sys.argv[1] == 'longest':
        orfs = longest_pep(pepfile)
        write_to_file(orfs)
    elif sys.argv[1] == 'clean':
        clean_defs(pepfile)
    else:
        print "longestORF.py [ longest | clean ]  [peptides.fa]"
