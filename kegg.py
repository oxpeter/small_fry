#!/usr/bin/env python
"""
A python module for dealing with KEGG pathways and orthologs
"""

import re
import sys

from itertools import chain
import argparse

from genomepy import config


class KeggTree(object):
    def __init__(self, kegg_f, cbir_ko):
        self.top_trees = {}                 # A level
        self.pathway_groups = {}            # B level
        self.pathway_groups_rev = {}        # B level
        self.pathway_groups_keggs_rev = {}  # B level
        self.pathways = {}                  # C level
        self.pathways_rev = {}              # C level

        self.pathways_terms = {}            # C level
        self.kegg_terms = {}                # D level

        self.cbir_dic = cbir_tree(cbir_ko) # for converting to cbir loci


        currentA = "none"
        aterm    = "none"
        currentB = "none"
        bterm    = "none"
        currentC = "none"
        cterm    = "none"
        currentD = "none"
        dterm    = "none"

        handle = open(kegg_f, 'rb')
        for line in handle:
            # regex returns (level), (term_id), (</b>), (definition)
            elements = re.search("^([ABCD])[\s<b>]+(\w+)(</b>)?\s?(.*)", line)
            if elements:
                if elements.group(1) == "A":
                    currentA = elements.group(2)

                    if currentA not in self.top_trees:
                        self.top_trees[currentA] = []

                elif elements.group(1) == "B":
                    currentB = elements.group(2)

                    if currentB not in self.pathway_groups:
                        self.pathway_groups[currentB] = []

                    self.top_trees[currentA] += [currentB]

                elif elements.group(1) == "C":
                    currentC = elements.group(2)
                    cterm = elements.group(4)
                    if currentC not in self.pathways:
                        self.pathways[currentC] = []
                        self.pathways_terms[currentC] = cterm

                    self.pathway_groups[currentB] += [currentC]

                    if currentC not in self.pathway_groups_rev:
                        self.pathway_groups_rev[currentC] = [currentB]
                    else:
                        self.pathway_groups_rev[currentC] += [currentB]

                elif elements.group(1) == "D":
                    currentD = elements.group(2)
                    dterm = elements.group(4)
                    self.kegg_terms[currentD] = dterm
                    self.pathways[currentC] += [currentD]

                    if currentD in self.pathways_rev:
                        self.pathways_rev[currentD] += [currentC]
                    else:
                        self.pathways_rev[currentD] = [currentC]

                    if currentD in self.pathway_groups_keggs_rev:
                        self.pathway_groups_keggs_rev[currentD] += [currentB]
                    else:
                        self.pathway_groups_keggs_rev[currentD] = [currentB]

    def convert_kegg(self, kegg):
        try:
            return self.cbir_dic[kegg]
        except KeyError:
            return None

    def find_pathways(self, kegg):
        "returns dictionary of all higher level pathways"
        paths = {'pathway_groups':[], 'pathways':[]}
        if kegg in self.pathway_groups_keggs_rev:
            paths['pathway_groups'] = self.pathway_groups_keggs_rev[kegg]
        if kegg in self.pathways_rev:
            paths['pathways'] = self.pathways_rev[kegg]
        if kegg in self.pathway_groups_rev: # not actually a kegg, but can still return
            paths['pathway_groups'] = self.pathway_groups_rev[kegg]
        return paths

    def list_keggs(self, pathway, terms=False, convert=False):
        if pathway in self.pathways:
            keggs = self.pathways[pathway]
        else:
            keggs = []
        if terms:
            keggs = [self.kegg_name(k) for k in keggs[:]]
        elif convert:
            keggs = list(chain.from_iterable(
                        [self.convert_kegg(k) for k in keggs[:] if self.convert_kegg(k)]
                        ))
        return keggs

    def list_pathways(self, pathway_group, terms=False):
        if pathway_group in self.pathway_groups:
            keggs = self.pathway_groups[pathway_group]
        else:
            keggs = []
        if terms:
            keggs = [self.pathway_name(k) for k in keggs[:]]

        return keggs

    def kegg_name(self, kegg):
        try:
            name = self.kegg_terms[kegg]
        except KeyError:
            name = "kegg not found"
        return name

    def pathway_name(self, pathway):
        try:
            name = self.pathways_terms[pathway]
        except KeyError:
            name = "pathway not found"
        return name

    def get_name(self, id):
        if id in self.pathways_terms:
            idname = self.pathways_terms[id]
            idtype = 'pathway'
        elif id in self.kegg_terms:
            idname = self.kegg_terms[kegg]
            idtype = 'kegg'
        return idname, idtype

    def search_terms(self, searchstring):
        results = {'keggs':[],'pathways':[],'pathway_groups':[],'top_level':[]}
        for kegg,term in self.kegg_terms.items():
            termsearch = re.search(searchstring,term)
            if termsearch:
                results['keggs'] += [kegg]
        for kegg,term in self.pathways_terms.items():
            termsearch = re.search(searchstring,term)
            if termsearch:
                results['pathways'] += [kegg]
        for kegg,term in self.pathway_groups.items():
            termsearch = re.search(searchstring,kegg)
            if termsearch:
                results['pathway_groups'] += [kegg]
        for kegg,term in self.top_trees.items():
            termsearch = re.search(searchstring,kegg)
            if termsearch:
                results['top_level'] += [kegg]
        return results

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "Performs various searches on the KEGG orthologs and pathways")
    # basic logging:
    parser.add_argument("-o", "--output", type=str, default='KEGG',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="don't print messages")

    # input options:
    parser.add_argument("-k", "--listkeggs", type=str,
                        help="""get KEGG orthologs for given KEGG pathway
                        Input KEGG path as number (e.g. 04612)""")
    parser.add_argument("-p", "--listpathways", type=str,
                        help="""list kegg pathways for given major pathway group
                        (e.g. Circulatory, Excretory, Immune, Transcription, Cell)""")
    parser.add_argument("--define", type=str,
                        help=""" Get the definition for a given list of KEGG orthologs.
                        """)
    parser.add_argument("-s", "--search", type=str,
                        help="search pathways for match to string")
    parser.add_argument("-l", "--level", type=str,
                        help="level for pathway search (A, B or C)")

    parser.add_argument("kegg_hierarchy", nargs=1,
                        help="""the KEGG hierarchy file.""")
    parser.add_argument("kegg_orthologs", nargs=1,
                        help="""the file showing the genes and their KEGG orthologs.
                        the file should contain two columns. The first is the gene name,
                        the second is the KEGG ortholog""")

    return parser

def cbir_tree(cbir_ko):
    handle = open(cbir_ko, 'rb')
    cbir_dic = {}
    for line in handle:
        cols = line.split()
        if len(cols) == 2:
            if cols[1] in cbir_dic:
                cbir_dic[cols[1]] += [cols[0]]
            else:
                cbir_dic[cols[1]] = [cols[0]]
    return cbir_dic

def load_config(config_f):
    config_h = open(config_f, 'rb')
    config_d = { line.split()[0]:line.split()[1] for line in config_h if line[0]!='#'}
    config_h.close()
    return config_d



if __name__ == '__main__':

    parser = define_arguments()
    args = parser.parse_args()

    kegg_f = args.kegg_hierarchy[0]   # location of .keg file (KEGG hierarchy)
    cbir_ko = args.kegg_orthologs[0] # location of gene - kegg ortholog list

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    kegg_tree = KeggTree(kegg_f, cbir_ko)


    #print "top_trees\n", kegg_tree.top_trees.items()[:5]                 # A level
    #print "pathway_groups\n", kegg_tree.pathway_groups.items()[:5]            # B level
    #print "pathway_groups_rev\n", kegg_tree.pathway_groups_rev.items()[:5]       # B level
    #print "pathway_groups_keggs_rev\n", kegg_tree.pathway_groups_keggs_rev.items()[:5]  # B level
    #print "pathways\n", kegg_tree.pathways.items()[:5]               # C level
    #print "pathways_rev\n", kegg_tree.pathways_rev.items()[:5]             # C level
    #print "pathways_terms\n", kegg_tree.pathways_terms.items()[:5]            # C level
    #print "kegg_terms\n", kegg_tree.kegg_terms.items()[:5]              # D level
    #print "cbir_dic\n", kegg_tree.cbir_dic.items()[:5] # for converting to cbir loci

    if args.listkeggs:
        for pathway in args.listkeggs.split(","):
            # get all kegg orthologs from given pathway:
            genes = kegg_tree.list_keggs(args.listkeggs, convert=False)

            verbalise("M",
                        "%d genes in pathway %s (%s)" %
                            (   len(genes),
                                args.listkeggs,
                                kegg_tree.pathway_name(args.listkeggs)
                            )
                    )
            verbalise("G", " ".join(genes))
            verbalise("G",
            "\n".join([ "%-8s %s" % (ko, kegg_tree.kegg_name(ko)) for ko in genes])
                )
    if args.listpathways:
        # get all pathways from given major pathway group:
        pathways = kegg_tree.list_pathways(args.listpathways)

        # print with definition:
        verbalise("Y",
            "\n".join([ "%-8s %s" % (p, kegg_tree.pathway_name(p)) for p in pathways])
                )

    if args.search:
        print args.search.split(",")
        for searchterm in args.search.split(","):
            if searchterm == "":
                continue
            print searchterm
            results = kegg_tree.search_terms(searchterm)
            for level in results:
                verbalise("M", level)
                if level == "pathways":

                    verbalise("Y",
            "\n".join([ "%-8s %s" % (p, kegg_tree.pathway_name(p)) for p in results[level]])
                        )
                elif level == "keggs":
                    verbalise("G",
            "\n".join([ "%-8s %s" % (ko, kegg_tree.kegg_name(ko)) for ko in results[level]])
                            )
                else:
                    verbalise("C", "\n".join(results[level]))

    if args.define:
        kos = args.define.split(',')
        for ko in kos:
            if ko == "":
                continue

            ko_def = kegg_tree.kegg_name(ko)
            if ko in kegg_tree.pathways_rev:
                higher_paths = kegg_tree.pathways_rev[ko]
            else:
                higher_paths = []
            hp_defs =  [ kegg_tree.pathway_name(p) for p in higher_paths ]
            hp_zip = zip(higher_paths, hp_defs)

            verbalise("G",  "%-8s %s" % (ko, ko_def) )
            verbalise("Y",  "\n".join(["%-8s %s" % (p, d) for p,d in hp_zip]))


