#!/usr/bin/env python

import argparse
import os
import re
import sys
import tempfile

import matplotlib.pyplot as plt

from ortholotree import config

############################################
def define_arguments():
    parser = argparse.ArgumentParser(description=
            "Creates pie charts of the number of taxa represented in a genome.")

    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='taxa_pie_chart',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")

    # data file options:
    parser.add_argument("blast_results", type=str, nargs='+',
                        help="all blast output files to combine and report on")
    parser.add_argument("--swissprot_definitions", type=str,
                        default= '/Volumes/antqueen/booster/PRO_Odontomachus/trinity_denovo_normalized_odontomachus/speclist.txt',
                        help="""location of swissprot taxa definition file.""")
    parser.add_argument("--taxa_groups", type=str,
                        default='/Volumes/antqueen/booster/PRO_Odontomachus/swissprot_classes.csv',
                        help="""location of the swissprot taxa classification file.""")

    # select ant subfamily:
    parser.add_argument("--ponerines", action='store_true',
                        help="""subpartition ponerine ants""")
    parser.add_argument("--myrmecines", action='store_true',
                        help="""subpartition myrmecine ants""")
    parser.add_argument("--formicines", action='store_true',
                        help="""subpartition formicine ants""")
    parser.add_argument("--ants", action='store_true',
                        help="""do not subpartition ants""")

    return parser

# test files used in original program:
"TD_predicted.blast_hymenoptera.out"
"blastp.outfmt6"


def get_taxa_groups():
    # define classes:
    taxa_groups = {
        "AGEAP":"Arthropods",
        "AGEOR":"Arthropods",
        "ARTSA":"Arthropods",
        "ARTSF":"Arthropods",
        "ARTSX":"Arthropods",
        "CARRO":"Arthropods",
        "CTEFE":"Arthropods",
        "CUPSA":"Arthropods",
        "HOMAM":"Arthropods",
        "LIMPO":"Arthropods",
        "PANAR":"Arthropods",
        "TACTR":"Arthropods",
        "PIG":"Bilateria",
        "RAT":"Bilateria",
        "AILME":"Bilateria",
        "ANAPL":"Bilateria",
        "ANSAN":"Bilateria",
        "ASTAS":"Bilateria",
        "BOVIN":"Bilateria",
        "BRAFL":"Bilateria",
        "CALJA":"Bilateria",
        "CANLF":"Bilateria",
        "CAVPO":"Bilateria",
        "CHICK":"Bilateria",
        "CHLAE":"Bilateria",
        "CRIGR":"Bilateria",
        "CYPCA":"Bilateria",
        "DANRE":"Bilateria",
        "DASNO":"Bilateria",
        "DIDAL":"Bilateria",
        "DIPOM":"Bilateria",
        "ELEEL":"Bilateria",
        "FELCA":"Bilateria",
        "GALSE":"Bilateria",
        "GORGO":"Bilateria",
        "HALRO":"Bilateria",
        "HELSU":"Bilateria",
        "HORSE":"Bilateria",
        "HUMAN":"Bilateria",
        "LITCT":"Bilateria",
        "LOPSP":"Bilateria",
        "MACFA":"Bilateria",
        "MACMU":"Bilateria",
        "MACNE":"Bilateria",
        "MESAU":"Bilateria",
        "MOUSE":"Bilateria",
        "MUSMC":"Bilateria",
        "MUSPF":"Bilateria",
        "MUSSP":"Bilateria",
        "MUSTR":"Bilateria",
        "MYOCO":"Bilateria",
        "NYCCO":"Bilateria",
        "ONCMY":"Bilateria",
        "OPSTA":"Bilateria",
        "PANTR":"Bilateria",
        "PAPHA":"Bilateria",
        "PAROL":"Bilateria",
        "PONAB":"Bilateria",
        "PONPY":"Bilateria",
        "PSETT":"Bilateria",
        "RABIT":"Bilateria",
        "RHIMB":"Bilateria",
        "SALSA":"Bilateria",
        "SHEEP":"Bilateria",
        "SPECI":"Bilateria",
        "SQUAC":"Bilateria",
        "TAEGU":"Bilateria",
        "TAKRU":"Bilateria",
        "THANI":"Bilateria",
        "TREBE":"Bilateria",
        "TRIHK":"Bilateria",
        "XENLA":"Bilateria",
        "XENTR":"Bilateria",
        "ASCSU":"Ecdysozoa",
        "CAEBR":"Ecdysozoa",
        "CAEEL":"Ecdysozoa",
        "ONCVO":"Ecdysozoa",
        "OSCTI":"Ecdysozoa",
        "RSVP":"Eukaryotes",
        "ACHKL":"Eukaryotes",
        "ACRMI":"Eukaryotes",
        "AGABT":"Eukaryotes",
        "ARATH":"Eukaryotes",
        "ASHGO":"Eukaryotes",
        "ASPFU":"Eukaryotes",
        "BLAAD":"Eukaryotes",
        "BOTFU":"Eukaryotes",
        "CANAL":"Eukaryotes",
        "CRYCU":"Eukaryotes",
        "CRYNJ":"Eukaryotes",
        "DICDI":"Eukaryotes",
        "EMENI":"Eukaryotes",
        "GIAIN":"Eukaryotes",
        "LEIIN":"Eukaryotes",
        "LEIMA":"Eukaryotes",
        "MAIZE":"Eukaryotes",
        "MEDSA":"Eukaryotes",
        "METMP":"Eukaryotes",
        "NAEGR":"Eukaryotes",
        "NEUCR":"Eukaryotes",
        "NEUIN":"Eukaryotes",
        "OENBE":"Eukaryotes",
        "ORYSJ":"Eukaryotes",
        "PARTE":"Eukaryotes",
        "PHAAO":"Eukaryotes",
        "PHYPO":"Eukaryotes",
        "PICGU":"Eukaryotes",
        "PLAF7":"Eukaryotes",
        "PLAYO":"Eukaryotes",
        "PLESA":"Eukaryotes",
        "PODAS":"Eukaryotes",
        "RHISN":"Eukaryotes",
        "SCHPO":"Eukaryotes",
        "SOLTU":"Eukaryotes",
        "TETTH":"Eukaryotes",
        "THACU":"Eukaryotes",
        "THETK":"Eukaryotes",
        "TOBAC":"Eukaryotes",
        "USTMA":"Eukaryotes",
        "YEAST":"Eukaryotes",
        "ZYMMO":"Eukaryotes",
        "Ador":"Hymenoptera",
        "Aflo":"Hymenoptera",
        "Amel":"Hymenoptera",
        "Aros":"Hymenoptera",
        "Bimp":"Hymenoptera",
        "Bter":"Hymenoptera",
        "Ccin":"Hymenoptera",
        "CoFl":"Hymenoptera",
        "Csol":"Hymenoptera",
        "Dall":"Hymenoptera",
        "Fari":"Hymenoptera",
        "Mdem":"Hymenoptera",
        "Mrot":"Hymenoptera",
        "Nlec":"Hymenoptera",
        "Nvit":"Hymenoptera",
        "Oabi":"Hymenoptera",
        "Pcan":"Hymenoptera",
        "Pdom":"Hymenoptera",
        "Tpre":"Hymenoptera",
        "APIME":"Hymenoptera",
        "LYSTE":"Hymenoptera",
        "NASVI":"Hymenoptera",
        "AEDAE":"Insecta",
        "ANOAL":"Insecta",
        "ANOGA":"Insecta",
        "ANOQU":"Insecta",
        "ANOST":"Insecta",
        "APILI":"Insecta",
        "BACDO":"Insecta",
        "BLACR":"Insecta",
        "BLADI":"Insecta",
        "BLAGE":"Insecta",
        "BOMIG":"Insecta",
        "BOMMO":"Insecta",
        "BOMPE":"Insecta",
        "BRACO":"Insecta",
        "CAMAT":"Insecta",
        "CARGR":"Insecta",
        "CERCA":"Insecta",
        "CERCO":"Insecta",
        "CHITE":"Insecta",
        "CHOFU":"Insecta",
        "CICVR":"Insecta",
        "CLOAL":"Insecta",
        "CULQU":"Insecta",
        "DROAN":"Insecta",
        "DROER":"Insecta",
        "DROEZ":"Insecta",
        "DROFU":"Insecta",
        "DROGR":"Insecta",
        "DROGU":"Insecta",
        "DROMA":"Insecta",
        "DROME":"Insecta",
        "DROMI":"Insecta",
        "DROMO":"Insecta",
        "DROPE":"Insecta",
        "DROPS":"Insecta",
        "DROSE":"Insecta",
        "DROSI":"Insecta",
        "DROTE":"Insecta",
        "DROVI":"Insecta",
        "DROWI":"Insecta",
        "DROYA":"Insecta",
        "GLOMM":"Insecta",
        "HELAM":"Insecta",
        "HELVI":"Insecta",
        "HOLDI":"Insecta",
        "LOCMI":"Insecta",
        "LONON":"Insecta",
        "LUCCR":"Insecta",
        "LUCCU":"Insecta",
        "LUCMI":"Insecta",
        "MANSE":"Insecta",
        "MUSDO":"Insecta",
        "OCHTR":"Insecta",
        "PERAM":"Insecta",
        "PLOIN":"Insecta",
        "PLUXY":"Insecta",
        "POLVA":"Insecta",
        "POPJA":"Insecta",
        "PROTE":"Insecta",
        "SCHAM":"Insecta",
        "SCHGA":"Insecta",
        "SCHGR":"Insecta",
        "SOLNG":"Insecta",
        "SOLS1":"Insecta",
        "SPOFR":"Insecta",
        "SPOLI":"Insecta",
        "TENMO":"Insecta",
        "TIMBA":"Insecta",
        "TRICA":"Insecta",
        "TRINI":"Insecta",
        "VESCR":"Insecta",
        "FFV":"Life",
        "BDVV":"Life",
        "KORV":"Life",
        "SFV1":"Life",
        "WDSV":"Life",
        "ACIAD":"Life",
        "ALKOO":"Life",
        "ANAMA":"Life",
        "AQUAE":"Life",
        "AVIRE":"Life",
        "AZOBR":"Life",
        "AZOC5":"Life",
        "AZOVI":"Life",
        "BACHD":"Life",
        "BACSU":"Life",
        "BPPHS":"Life",
        "BPS13":"Life",
        "BRUSU":"Life",
        "CALS8":"Life",
        "CAUCR":"Life",
        "CHRVO":"Life",
        "CILVC":"Life",
        "CLOK5":"Life",
        "COXBU":"Life",
        "CROS8":"Life",
        "CYLCV":"Life",
        "ECOLI":"Life",
        "ECOLX":"Life",
        "EHRCJ":"Life",
        "EHRCR":"Life",
        "FOAMV":"Life",
        "GLUOY":"Life",
        "GLVWB":"Life",
        "HAEIF":"Life",
        "HAEIN":"Life",
        "HELPY":"Life",
        "HERAU":"Life",
        "HV1B1":"Life",
        "HV1MP":"Life",
        "HV1ND":"Life",
        "HV2D1":"Life",
        "HV2D2":"Life",
        "HV2RO":"Life",
        "HV2SB":"Life",
        "HV2UC":"Life",
        "JEMBR":"Life",
        "LACS1":"Life",
        "LACSN":"Life",
        "LISMO":"Life",
        "MIMIV":"Life",
        "MLVCB":"Life",
        "MLVMS":"Life",
        "MLVRK":"Life",
        "NEIGO":"Life",
        "NORAV":"Life",
        "NOSS1":"Life",
        "NPVAC":"Life",
        "NPVLD":"Life",
        "NPVOP":"Life",
        "NPVSL":"Life",
        "PARL1":"Life",
        "PCV87":"Life",
        "PECAS":"Life",
        "PROMA":"Life",
        "PROMI":"Life",
        "PSYCK":"Life",
        "PVCV2":"Life",
        "RHIL3":"Life",
        "RHIME":"Life",
        "RHIS3":"Life",
        "RHOPA":"Life",
        "RHORH":"Life",
        "RICBR":"Life",
        "RICCN":"Life",
        "RICFE":"Life",
        "RICPR":"Life",
        "ROSS1":"Life",
        "SBWMN":"Life",
        "SFV3L":"Life",
        "SFVCP":"Life",
        "SHESR":"Life",
        "SIVG1":"Life",
        "SIVMB":"Life",
        "SIVVG":"Life",
        "STAES":"Life",
        "STAMF":"Life",
        "SYNY3":"Life",
        "WOLPI":"Life",
        "WOLPM":"Life",
        "WOLPP":"Life",
        "WOLSP":"Life",
        "WOLTR":"Life",
        "WOLWR":"Life",
        "ANTEL":"Metazoa",
        "CIOIN":"Metazoa",
        "HELCR":"Metazoa",
        "HEMPU":"Metazoa",
        "HYDVU":"Metazoa",
        "LYTPI":"Metazoa",
        "MARGL":"Metazoa",
        "NEMVE":"Metazoa",
        "PATPE":"Metazoa",
        "STRIE":"Metazoa",
        "STRPU":"Metazoa",
        "TIGCA":"Metazoa",
        "TRIGR":"Metazoa",
        "APLCA":"Protostome",
        "CRAGI":"Protostome",
        "LYMST":"Protostome",
        "SABMA":"Protostome",
        "SCHHA":"Protostome",
        "SCHMA":"Protostome",
        }
    if args.ponerines:
        taxa_groups.update({ "Dqua":"Ponerines", "Hsal":"Ponerines",
                            "Acep":"Ants", "Aech":"Ants", "Cbir":"Ants", "Cflo":"Ants",
                            "Ebur":"Ants", "Fexs":"Ants", "Lhum":"Ants",
                            "Mpha":"Ants", "Pbar":"Ants", "Sinv":"Ants", "Veme":"Ants",
                            "Waur":"Ants", "CAMFO":"Ants", "SOLIN":"Ants",
                             })

    elif args.myrmecines:
        taxa_groups.update({ "Acep":"Myrmecines", "Aech":"Myrmecines", "Cbir":"Ants", "Cflo":"Ants",
                            "Ebur":"Ants", "Fexs":"Ants", "Hsal":"Ants", "Lhum":"Ants",
                            "Mpha":"Myrmecines", "Pbar":"Myrmecines", "Sinv":"Myrmecines",
                            "Veme":"Myrmecines", "Dqua":"Ants",
                            "Waur":"Myrmecines", "CAMFO":"Ants", "SOLIN":"Myrmecines",
                             })
    elif args.formicines:
        taxa_groups.update({ "CAMFO":"Formicines", "Cflo":"Formicines", "Fexs":"Formicines",
                            "Acep":"Ants", "Aech":"Ants", "Cbir":"Ants",
                            "Ebur":"Ants", "Hsal":"Ants", "Lhum":"Ants",
                            "Mpha":"Ants", "Pbar":"Ants", "Sinv":"Ants", "Veme":"Ants",
                            "Waur":"Ants", "SOLIN":"Ants", "Dqua":"Ants",
                             })

    elif args.ants:
        taxa_groups.update({ "Acep":"Ants", "Aech":"Ants", "Cbir":"Ants", "Cflo":"Ants",
                            "Ebur":"Ants", "Fexs":"Ants", "Hsal":"Ants", "Lhum":"Ants",
                            "Mpha":"Ants", "Pbar":"Ants", "Sinv":"Ants", "Veme":"Ants",
                            "Waur":"Ants", "CAMFO":"Ants", "SOLIN":"Ants", "Dqua":"Ants",
                             })
    else: # this option is currently the same as args.ants, and is added just to have a default
        taxa_groups.update({ "Acep":"Ants", "Aech":"Ants", "Cbir":"Ants", "Cflo":"Ants",
                            "Ebur":"Ants", "Fexs":"Ants", "Hsal":"Ants", "Lhum":"Ants",
                            "Mpha":"Ants", "Pbar":"Ants", "Sinv":"Ants", "Veme":"Ants",
                            "Waur":"Ants", "CAMFO":"Ants", "SOLIN":"Ants", "Dqua":"Ants",
                             })

    return taxa_groups




def make_autopct(values):
    """
    This function is used by plt.pie to return only values greater than 1%, to reduce
    clutter on the graph.
    """
    def my_autopct(pct):
        if pct < 1:
            return ''
        else:
            return "%.1f%%" % pct
    return my_autopct

def pie_chart(taxcount, filename="pie_chart.pdf", truncate=True):

    # collate all small classes (< 1%) together if truncate == True
    total_trans = sum(taxcount.values())
    if truncate:
        orig_keys = taxcount.keys()[:]

        # add new class to combine old results into:
        newcount = { 'others':0 }
        for t in orig_keys:
            if float(taxcount[t]) / total_trans < 0.01:
                newcount['others'] += taxcount[t]
            else:
                newcount[t] = taxcount[t]

        # because of collation, we can print all % values
        autopct = "%.1f%%"

    else:
        newcount = taxcount
        autopct = make_autopct(newcount.values())


    plt.figure(figsize=(10,5))
    #ax1 = plt.subplot2grid( (1,2),(0,0))
    ax1 = plt.subplot(111)
    #ax1.set_axis_bgcolor('white')



    pie_wedge_collection = ax1.pie(newcount.values(),
                                    labels=newcount.keys(),
                                    labeldistance=1.05,
                                    autopct=autopct)

    for pie_wedge in pie_wedge_collection[0]:
        pie_wedge.set_edgecolor('white')

    # set axis ratios to make pie a circle
    plt.axis('equal')
    plt.title("Taxa distribution of best blast hits\n(%d transcripts)" % (total_trans))
    #ax2 = plt.subplot2grid( (10,5),(0,5), rowspan=5, colspan=5)

    ax1.legend(loc=3)
    plt.tight_layout()
    plt.savefig(filename, format='pdf')
    plt.show()

if __name__ == '__main__':
    dbpaths = config.import_paths()

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    temp_dir = tempfile.mkdtemp()
    os.rmdir(temp_dir)  # dir must be empty!

    best_hit = {}

    """
    # populate with hymenopteran hits
    handle = open(hymen, 'rb')
    for line in handle:
        gene = line.split()[0]
        score = float(line.split()[2])
        species = line.split()[1]
        best_hit[gene] = (score, species)
    """


    # add all genes to the list, replacing existing entries if the score is higher.
    for bf in args.blast_results:
        handle = open(bf, 'rb')
        for line in handle:
            score = float(line.split()[2])
            gene = line.split()[0]
            species = line.split()[1]
            if gene in best_hit:
                if score > best_hit[gene][0]:
                    best_hit[gene] = (score, species)
            else:
                best_hit[gene] = (score, species)


    # create dictionary of swissprot species descriptors
    definitions = {}
    handle = open(args.swissprot_definitions, 'rb')

    for line in handle:
        if line.strip() !=  "" and line[0] != " ":
            definitions[line.split()[0]] = line.strip()
            prior = line.split()[0]
        else:
            definitions[prior] += line.strip()

    handle.close()

    # cluster by species
    results = {"unknown":0}

    for g in best_hit:
        if best_hit[g][0] < 70:
            results["unknown"] += 1
        else:
            metaz = best_hit[g][1].split('_')[-1]
            if metaz in definitions:
                if metaz in results:
                    results[metaz] += 1
                else:
                    results[metaz] = 1
            else:
                hymen = best_hit[g][1].split('|')[0]
                if hymen in results:
                    results[hymen] += 1
                else:
                    results[hymen] = 1

    for sp in results:
        if sp in definitions:
            print "%-6s\t%d\t%s" % (sp, results[sp], definitions[sp])
        else:
            print "%-6s\t%d" % (sp, results[sp])

    print '#' * 45

    # cluster species by taxonomic level (each being monophyletic with ponerine ants)
    taxconv = get_taxa_groups()

    taxcount = { t:0 for t in taxconv.values() }
    taxcount['unknown'] = 0

    for sp in results:
        if sp in taxconv:
            taxcount[taxconv[sp]] += results[sp]
        else:
            taxcount['unknown'] += results[sp]

    for t in taxcount:
        print "%-20s %d" % (t, taxcount[t])

    # plot all transcripts
    pie_chart(taxcount, logfile[:-3] + 'chart.pdf')

    # remove poor matches and plot again
    if 'unknown' in taxcount:
        del taxcount['unknown']
    pie_chart(taxcount, logfile[:-3] + 'accepted_chart.pdf' )


    ########################################################


