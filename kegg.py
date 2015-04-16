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
        if terms:
            keggs = [self.kegg_name(k) for k in keggs[:]]
        elif convert:
            keggs = list(chain.from_iterable(
                        [self.convert_kegg(k) for k in keggs[:] if self.convert_kegg(k)]
                        ))
        return keggs

    def list_pathways(self, pathway_group, terms=False):
        if pathway in self.pathway_groups:
            keggs = self.pathway_groups[pathway]
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

    # input options
    parser.add_argument("-k", "--listkeggs", type=str,
                        help="get Cbir genes for given KEGG pathway")
    parser.add_argument("-p", "--listpathways", type=str,
                        help="list kegg pathways for given pathway group")

    parser.add_argument("-s", "--search", type=str,
                        help="search pathways for match to string")
    parser.add_argument("-l", "--level", type=str,
                        help="level for pathway search (A, B or C)")

    parser.add_argument("-o", "--output", type=str, default='lorf2.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="don't print messages")

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
    if len(sys.argv)<2:
        print "Please supply the config file in the command line"
        exit()

    config_d = load_config(sys.argv[1])

    kegg_f = config_d['kegg_f']   # location of .keg file (KEGG hierarchy)
    cbir_ko = config_d['cbir_ko'] # location of gene - kegg ortholog list

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    kegg_tree = KeggTree(kegg_f, cbir_ko)
    #cbir_dic = cbir_tree(cbir_ko)

    if args.listkeggs:
        # get all keggs from given pathway:
        genes = kegg_tree.list_keggs(args.listkeggs, convert=True)
        print genes
        # convert to cbir genes:
        #genes = []
        #for k in keggs:
        #    try:
        #        genes += cbir_dic[k]
        #    except KeyError:
        #        pass
        verbalise("M", "%d genes in pathway %s (%s)" % ( len(genes), args.listkeggs, kegg_tree.pathway_name(args.listkeggs)))
        verbalise("G", "\n".join(genes))

    if args.listpathways:
        # get all pathways from given pathway group:
        pathways = kegg_tree.list_pathways(args.listpathways)
        verbalise("Y", "\n".join([ "%-8s %s" % (p, kegg_tree.pathway_name(p)) for p in pathways]))

    if args.search:
        results = kegg_tree.search_terms(args.search)
        for level in results:
            verbalise("M", level)
            verbalise("C", "\n".join(results[level]))

