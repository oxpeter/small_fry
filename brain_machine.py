#!/usr/bin/env python
"""
A preliminary script for working through machine learning of brain gene expression
data
"""

import re
import os
import sys

import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn import svm
from sklearn.feature_selection import RFE, RFECV
from sklearn import preprocessing

from genomepy import config
import kegg


verbalise = config.check_verbose()

###########  FUNCTIONS  ###############
def factorial(n):
    if n < 2: return 1
    return reduce(lambda x, y: x*y, xrange(2, int(n)+1))

def prob(s, p, n):
    " returns the cumulative probability P(X >= x) "
    s = int(s - 1)
    x = int(1.0 - p)
    a = n - s
    b = s + 1
    c = int(a + b - 1)
    prob = 0.0
    for j in xrange(a, c + 1):
        prob += factorial(c) / (factorial(j)*factorial(c-j)) \
                * x**j * (1 - x)**(c-j)
    return prob


def ncbi_dic(file):
	handle = open(file, 'rb')
	dic = { line.split()[0] : " ".join(line.split()[1:]) for line in handle }
	return dic

def show_name(loc, dic=None):
	if not dic:
		dic = ncbi_dic()
	try:
		return "%-14s %s" % (loc, dic[loc])
	except KeyError:
		return "%-14s --name not found--" % (loc)

def assemble_df(htseq_list, lognorm=True, sfnorm=False):
    "import and assemble raw counts into dataframes"
    dfs = {}
    for i,filename in enumerate(htseq_list):
        headname = os.path.basename(filename)[:-30]
        dfs[i] = pd.read_csv(filename, sep="\t", header=None, skipfooter=5, engine='python', index_col=0, names=[headname])
    df = pd.concat(dfs.values(), axis=1)
    if lognorm:
        # take the log2(x + 1) of each gene's expression (x)
        df=np.log2(df+1)
    if sfnorm:
        # divide each gene's expression value by the mean of all genes' values
        df=df/df.mean(axis=0)
    return df

def index_orthos(df, orthofile):
    orthodf = pd.read_csv(orthofile, sep="\t", header=None, index_col=0, names=['orthos'])
    newdf = df.join(orthodf).dropna()
    newdf.set_index(['orthos'], inplace=True)
    return newdf

def to_binary(names):
    binary = [ int(re.search("(queen|worker|foundress|gyne|_FL|_SP|Emerged)",name).group(1) in ['queen','foundress','gyne','_SP'] ) for name in names  ]
    return binary

def create_xy(df1, df2, normalise=False):
    """
    given two dataframes, find common loci, and return two numpy arrays for the
    genotypes and their corresponding y series of phenotypes (0 for worker/foraging, 1
    for queen/statary
    """
    jointdf = df1.join(df2).dropna()
    # remove duplicates:
    jointdf["index"] = jointdf.index
    jointdf.drop_duplicates(subset='index', take_last=True, inplace=True)
    del jointdf["index"]

    y1 = to_binary(df1.columns)
    y2 = to_binary(df2.columns)
    if normalise:
        scaled = preprocessing.scale(jointdf.T.values)
        X1 = scaled[:len(df1.columns), :]
        X2 = scaled[len(df1.columns):, :]
    else:
        X1 = jointdf.ix[:,:len(df1.columns)].T.values
        X2 = jointdf.ix[:,len(df1.columns):].T.values

    return X1, y1, X2, y2, jointdf.index

def predictor_svm(X1, y1, X2, y2):
    verbalise("C", "Size of training set: %d\nSize of test set: %d" % (len(y1), len(y2)))
    verbalise("C", "Number of features: %d" % X1.shape[1])
    clf = svm.SVC(kernel='linear')
    #print clf
    clf.fit(X1, y1)
    verbalise("Y", "real: %s\npred: %s" % (
            " ".join([str(x) for x in y2]), " ".join([ str(x) for x in clf.predict(X2)])
            ))
    cscore = clf.score(X2, y2)
    pval = prob( round(cscore*len(y2)), 0.5, len(y2) )
    verbalise("G", "score: %.2f (p=%.5f)\n" % (cscore, pval ))
    return clf.coef_, clf.score(X2, y2)

def get_features_rfe(X1, y1, X2, ffn=100):
    verbalise("C", "Sample size of set: %d" % (len(y1)) )
    verbalise("C", "Initial number of features: %d" % X1.shape[1] )
    verbalise("C", "Final number of features: %d" % ffn )
    estimator = svm.SVR(kernel="linear")
    selector = RFE(estimator, ffn, step=1)
    selector = selector.fit(X1, y1)
    newX1 = selector.transform(X1)
    newX2 = selector.transform(X2)
    return newX1, newX2, selector.support_

def get_features_rfecv(X1, y1, X2, step=1, cv=None):
    if cv:
        pass
    else:
        cv = len(y1)

    verbalise("C", "Sample size of set: %d" % (len(y1)) )
    verbalise("C", "Initial number of features: %d" % X1.shape[1] )
    estimator = svm.SVR(kernel="linear")
    selector = RFECV(estimator, cv=cv, step=step, scoring=None)
    selector = selector.fit(X1, y1)
    newX1 = selector.transform(X1)
    newX2 = selector.transform(X2)
    verbalise("C",
        "Final number of features: %d (average of %d CV scores = %.3f)" %
            (newX1.shape[1], cv, selector.grid_scores_[newX1.shape[1] - 1] ))
    return newX1, newX2, selector.support_, selector.grid_scores_

def multitest(df1, df2, name1='X1', name2='X2', ffn=100, normtogether=False, normind=False, cv=None):
    "normtogether normalises training and testing data together. "
    verbalise("M", "Predicting %s phase with %s using linear kernel SVM" % (name2, name1))
    X1, y1, X2, y2, dfindex = create_xy(df1, df2, normalise=normtogether)
    if normind:
        X1 = preprocessing.scale(X1)
        X2 = preprocessing.scale(X2)
    weights, df2_score = predictor_svm(X1, y1, X2, y2)
    verbalise("M", "Using RFE %s" % (name1))
    X1s, X2s, rfemask = get_features_rfe(X1, y1, X2, ffn=ffn)
    df2_score = predictor_svm(X1s, y1, X2s, y2)
    verbalise("M", "Using RFECV on %s" % (name1))
    X1c, X2c, rfecvmask, rfecv_gscores = get_features_rfecv(X1, y1, X2, cv=cv)
    df2_score = predictor_svm(X1c, y1, X2c, y2)
    print "-" * 75
    rfe_features = [ pair[0] for pair in zip(dfindex, rfemask) if pair[1] ]
    rfecv_features = [ pair[0] for pair in zip(dfindex, rfecvmask) if pair[1] ]

    return weights, rfe_features, rfecv_features, dfindex

def report_genes(weights, rfes, rfecvs, df_index, num=25):
    weighting = sorted(zip(df_index, weights[0]), key=lambda pair: pair[1], reverse=True)
    print "\n%d most weighted features from SVM (positive and negative): " % num
    verbalise("G", '\n'.join([ "%+7.4f %s" % (tup[1],show_name(tup[0],ncbi))
                                for tup in weighting[:num]]))
    verbalise("C", '\n'.join([ "%+7.4f %s" % (tup[1],show_name(tup[0],ncbi))
                                for tup in weighting[-num:]]))
    print "\n%d features from recursive feature elimination: (limit set by user)" % len(rfes)
    verbalise("Y", '\n'.join([ show_name(loc, ncbi) for loc in rfes]))
    print "\n%d features from recursive feature elimination with cross validation: " % len(rfecvs)
    verbalise("B", '\n'.join([ show_name(loc, ncbi) for loc in rfecvs]))

def load_config(config_f):
    config_h = open(config_f, 'rb')
    config_d = { line.split()[0]:line.split()[1] for line in config_h if line[0]!='#'}
    config_h.close()
    return config_d


###########  ANALYSES  ################
if __name__ == '__main__':
    if len(sys.argv)<2:
        print "Please supply the config file in the command line"
        exit()

    config_d = load_config(sys.argv[1])

    # create lists of all DEGs from each experiment:
    meth_degs = config.make_a_list(config_d['meth_degs_f'])
    brsw_degs = config.make_a_list(config_d['brsw_degs_f'])
    acro_degs = config.make_a_list(config_d['acro_degs_f'])

    # create dictionary for extracting ncbi names of Cbir genes:
    ncbi = ncbi_dic(file=config_d['ncbi_names'])

    # collect all raw count data:
    meth_r = [ os.path.join(config_d['meth_rd'],file)
                for file in os.listdir(config_d['meth_rd']) if re.search("^M.*htseq.gene",file) ]
    brsw_r   = [ os.path.join(config_d['brsw_rd'],file)
                for file in os.listdir(config_d['brsw_rd']) if re.search("^[RP].*htseq.gene",file) ]
    polk_r = [ os.path.join(config_d['polk_rd'],file)
                for file in os.listdir(config_d['polk_rd']) if re.search("htseq.gene",file) ]
    acro_r = [ os.path.join(config_d['acro_rd'],file)
                for file in os.listdir(config_d['acro_rd']) if re.search("htseq.gene",file) ]
    polu_r = [ os.path.join(config_d['polu_rd'],file)
                for file in os.listdir(config_d['polu_rd']) if re.search("htseq.gene",file) ]
    bsse_r = [ os.path.join(config_d['bsse_rd'],file)
                for file in os.listdir(config_d['bsse_rd']) if re.search("(_FL|_SP).+htseq.gene",file) ]

    #verbalise( "Number of files found for each study:")
    #verbalise("G", "Meth: %d\nBrSw: %d\nBSse: %d\nPolK: %d\nPolU: %d\nAcro: %d" % (len(meth_r), len(brsw_r), len(bsse_r), len(polk_r), len(polu_r), len(acro_r)))

    print "\nAssembling dataframes...\n"

    meth_df = assemble_df(meth_r, lognorm=True, sfnorm=False)
    brsw_df = assemble_df(brsw_r, lognorm=True, sfnorm=False)
    bsse_df = assemble_df(bsse_r, lognorm=True, sfnorm=False)
    polk_df = assemble_df(polk_r, lognorm=True, sfnorm=False)
    polu_df = assemble_df(polu_r, lognorm=True, sfnorm=False)
    acro_df = assemble_df(acro_r, lognorm=True, sfnorm=False)
    cbir_df = meth_df.join(brsw_df)#.join(bsse_df)

    # create dataframes of vst:
    meth_vst = pd.read_csv(config_d['meth_f'], sep=" ", header=0)
    brsw_vst = pd.read_csv(config_d['brsw_f'], sep=" ", header=0)
    polu_vst = pd.read_csv(config_d['polu_f'], sep=" ", header=0)
    polk_vst = pd.read_csv(config_d['polk_f'], sep=" ", header=0)
    acro_vst = pd.read_csv(config_d['acro_f'], sep=" ", header=0)

    # convert to orthologs:
    polu_ortho = index_orthos(polu_df, config_d['convert_orthos'])
    polk_ortho = index_orthos(polk_df, config_d['convert_orthos'])
    acro_ortho = index_orthos(acro_df, config_d['convert_orthos'])

    all_degs = meth_degs.keys()+brsw_degs.keys()

    # prepare Kegg objects:
    kegg_tree = kegg.KeggTree(config_d['kegg_f'],config_d['cbir_ko'])



    ###########################  TESTING SITE  #############################

    #genelist = list(set(["LOC105281428_Q","LOC105278524","LOC105280170","LOC105282036","LOC105275236","LOC105277456","LOC105283266","LOC105284141","LOC105285273","LOC105281340","LOC105286142","LOC105284539","LOC105279365","LOC105279366","LOC105275821","LOC105283812","LOC105287615","LOC105277435","LOC105287013","LOC105282105","LOC105285370","LOC105288082","LOC105287763","LOC105281770","LOC105285920","LOC105281330","LOC105279396","LOC105279397","LOC105285921","LOC105284974","LOC105275269","LOC105285466","LOC105280465","LOC105283721","LOC105275039","LOC105277654","LOC105275039","LOC105277654","LOC105285466","LOC105285597","LOC105288201","LOC105281757","LOC105277039","LOC105284284","LOC105287193","LOC105280785"]))
    genelist = list(cbir_df.index[:])
    #genelist = list(set(brsw_degs) - (set(meth_degs) & set(brsw_degs)))
    #genelist = list((set(brsw_degs) | set(meth_degs)) - (set(meth_degs) & set(brsw_degs)))
    #genelist = list(set(brsw_degs))


    # Generate random datasets to test performance:
    genelist = list(cbir_df.index[:])

    all_scores = []
    for i in range(1000):
        random.shuffle(genelist)
        random_set = genelist[:5]
        X1, y1, X2, y2, dfindex = create_xy(cbir_df.loc[random_set], bsse_df)
        weight, score = predictor_svm(X1, y1, X2, y2)
        all_scores.append(score)

    plt.hist(all_scores, bins=24)
    plt.show()

    ###########################  WINNING ANALYSES  #########################

    """
    pairwisecombinations=[(brsw_df.loc[brsw_degs],meth_df,"Broodswap","Methylation PE",6),
    (brsw_df.loc[brsw_degs],bsse_df,"Broodswap","Broodswap SE",6),
    (brsw_df.loc[brsw_degs],acro_ortho,"Broodswap","Acromyrmex heads",6),
    (brsw_df.loc[brsw_degs],polk_ortho,"Broodswap","Polistes F13",6),
    (meth_df.loc[meth_degs],brsw_df,"Methylation","Broodswap",4),
    (meth_df.loc[meth_degs],bsse_df,"Methylation","Broodswap SE",4),
    (meth_df.loc[meth_degs],acro_ortho,"Methylation","Acromyrmex heads",4),
    (meth_df.loc[meth_degs],polk_ortho,"Methylation","Polistes F13",4),
    (cbir_df.loc[set(brsw_degs.keys() + meth_degs.keys())],bsse_df,"Broodswap + Methylation","Broodswap SE",5),
    (cbir_df.loc[set(brsw_degs.keys() + meth_degs.keys())],acro_ortho,"Broodswap + Methylation","Acromyrmex heads",5),
    (cbir_df.loc[set(brsw_degs.keys() + meth_degs.keys())],polk_ortho,"Broodswap + Methylation","Polistes F13",5),
    (index_orthos(acro_df.loc[acro_degs], convert_orthos), cbir_df.loc[set(brsw_degs.keys() + meth_degs.keys())], "Acromyrmex heads", "Cerapachys controls",3),
    ]

    for df1, df2, name1, name2, cv in pairwisecombinations:
        weights, rfes, rfecvs, dfindex = multitest(
                df1, df2,
                name1,name2,
                ffn=5, cv=cv, normtogether=False, normind=True)
        report_genes(weights, rfes, rfecvs, dfindex, 10)
    """


    """
    # Generate random datasets to test performance:
    genelist = list(cbir_df.index[:])

    all_scores = []
    for i in range(1000):
        random.shuffle(genelist)
        random_set = genelist[:5]
        X1, y1, X2, y2, dfindex = create_xy(cbir_df.loc[random_set], bsse_df)
        weight, score = predictor_svm(X1, y1, X2, y2)
        all_scores.append(score)

    plt.hist(all_scores, bins=24)
    plt.show()
    """


    """
    # do machine learning on KEGG pathway:
    pathway = '00562'
    genes = kegg_tree.list_keggs(pathway, convert=True)
    verbalise("Y", "#"*60, "\n%d genes found in %s" % (len(genes), kegg_tree.pathway_name(pathway)))
    weights, rfes, rfecvs, dfindex = multitest(
                                    cbir_df.loc[genes], bsse_df,
                                    kegg_tree.pathway_name(pathway),'broodswap SE',
                                    ffn=5, cv=5, normtogether=False, normind=True)
    """

    """
    # test all KEGG pathways:
    pathways = kegg_tree.pathways.keys()
    bestperformers = []
    for path in pathways:
        genes = kegg_tree.list_keggs(path, convert=True)
        genelist = list(set(genes))

        print "testing pathway %s:" % (kegg_tree.pathway_name(path))
        print "%d genes in initial genelist." % len(genelist)

        all_scores = []
        for i in range(1000):
            random.shuffle(genelist)
            random_set = genelist[:5]
            X1, y1, X2, y2, dfindex = create_xy(cbir_df.loc[random_set], bsse_df)
            #scaler = preprocessing.StandardScaler()
            #X1 = scaler.fit_transform(X1)
            #X2 = scaler.transform(X2)
            X1 = preprocessing.scale(X1)
            X2 = preprocessing.scale(X2)
            weight, score = predictor_svm(X1, y1, X2, y2)
            all_scores.append(score)
        print "mean:   %.2f\tmedian: %.2f\n" % (np.mean(all_scores), np.median(all_scores))
        if np.median(all_scores) > 0.75:
            bestperformers.append((path, kegg_tree.pathway_name(path)))

    for pathid, pathname in bestperformers:
        print pathid, pathname
    """


    ########################################################################


