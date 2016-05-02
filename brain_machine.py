#!/usr/bin/env python
"""
A preliminary script for working through machine learning of brain gene expression
data
"""

import re
import os
import sys
import argparse

import random
import pandas as pd
import numpy as np
from numpy.linalg import norm
from scipy.stats import entropy
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn import svm
from sklearn.feature_selection import RFE, RFECV
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier

from genomepy import config
import kegg





###########  FUNCTIONS  ###############
def define_arguments():
    parser = argparse.ArgumentParser(description=
            "performs various support vector machine learning algorithms on RNA-Seq datasets")

    # input options
    parser.add_argument("configfile", metavar='filename', nargs=1, type=str,
                        help="the config file with links to all the data files")

    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='lorf2.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")

    parser.add_argument("-M", "--multitest", action='store_true',default=False,
                        help="perform three SVM analyses on each set")
    parser.add_argument("-f", "--ffn", type=int, default=5,
                        help="final feature number for recursive feature elimination")
    parser.add_argument("-c", "--cv", type=int, default=10,
                        help="cross validation number")

    parser.add_argument("-R", "--randomsets", action='store_true',default=False,
                        help="iterate through random subsets of each set of genes")
    parser.add_argument("-r", "--replicates", type=int, default=1000,
                        help="number of replicates (default=1000)")
    parser.add_argument("-g", "--genenum", type=int, default=5,
                        help="number of genes to chose from complete list (default=5)")


    return parser

def JSD(P, Q):
    """
    calculates the Jenson Shannon Divergence. Code from Doug Shore
    http://stackoverflow.com/questions/15880133/jensen-shannon-divergence
    """
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return 0.5 * (entropy(_P, _M) + entropy(_Q, _M))

def factorial(n):
    if n < 2: return 1
    return reduce(lambda x, y: x*y, xrange(2, int(n)+1))

def prob(s, p, n):
    " returns the cumulative probability P(X >= x) "
    s = (s - 1)
    x = (1.0 - p)
    a = int(n - s)
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

def sort_orthos(row):
    if row['ortho1'].str.split('|').str[0] == 'Cbir':
        return row['ortho2'].str.split('|').str[1], row['ortho1'].str.split('|').str[1]
    elif row['ortho2'].str.split('|').str[0] == 'Cbir':
        return row['ortho1'].str.split('|').str[1], row['ortho2'].str.split('|').str[1]
    else:
        return np.nan, np.nan

def create_ortho_df(orthofile):
    # create new dataframe from orthomcl ortholog output file:
    orthodf = pd.read_csv(orthofile, sep="\t", header=None, names=['ortho1','ortho2', 'score'])

    # create new column for Cerapachys orthologs:
    orthodf['cbir_loc'] = np.nan
    orthodf.loc[ orthodf['ortho1'].str.split('|').str[0] == 'Cbir', 'cbir_loc' ] = orthodf['ortho1'].str.split('|').str[1]
    orthodf.loc[ orthodf['ortho2'].str.split('|').str[0] == 'Cbir', 'cbir_loc' ] = orthodf['ortho2'].str.split('|').str[1]
    # create new column for other species' orthologs:
    orthodf = orthodf.dropna(subset = ['cbir_loc'])
    orthodf['other_loc'] = np.nan
    orthodf.loc[ orthodf['ortho1'].str.split('|').str[0] == 'Cbir', 'other_loc' ] = orthodf['ortho2'].str.split('|').str[1]
    orthodf.loc[ orthodf['ortho2'].str.split('|').str[0] == 'Cbir', 'other_loc' ] = orthodf['ortho1'].str.split('|').str[1]

    # set index of df to the non-Cbir ortholog:
    orthodf.set_index(['other_loc'], inplace=True)
    orthodf = orthodf[['cbir_loc']]
    return orthodf

def index_orthos(df, orthodf):
    # join to other dataframe, which is indexed with non-Cbir orthologs:
    newdf = df.join(orthodf).dropna()  # also removes any rows with no orthologs
    newdf.set_index(['cbir_loc'], inplace=True)
    return newdf

def to_binary(names):
    try:
        binary = [ int(re.search("(QUEEN|WORKER|FOUNDRESS|GYNE|FORAGER|NURSE|_FL|_SP|EMERGED|ALPHA|LOW)",name.upper()).group(1) in ['QUEEN','FOUNDRESS','GYNE','_SP','ALPHA'] ) for name in names  ]
    except AttributeError:
        print names
        exit()

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
    pval = prob( len(y2)-round(cscore*len(y2)) + 1, 0.5, len(y2) )
    verbalise("G", "score: %.2f (p=%.5f)\n" % (cscore, pval ))
    return clf.coef_, clf.score(X2, y2)

def random_forest_engine(X1, y1, X2, y2):
    verbalise("C", "Size of training set: %d\nSize of test set: %d" % (len(y1), len(y2)))
    verbalise("C", "Number of features: %d" % X1.shape[1])
    clf = RandomForestClassifier(max_features='log2', n_estimators=120, random_state=125)
    #print clf
    clf.fit(X1, y1)
    verbalise("Y", "real: %s\npred: %s" % (
            " ".join([str(x) for x in y2]), " ".join([ str(x) for x in clf.predict(X2)])
            ))
    cscore = clf.score(X2, y2)
    pval = prob( len(y2)-round(cscore*len(y2)) + 1, 0.5, len(y2) )
    verbalise("G", "score: %.2f (p=%.5f)\n" % (cscore, pval ))
    return clf.score(X2, y2)

def get_features_rfe(X1, y1, X2, ffn=100):
    verbalise("C", "Sample size of set: %d" % (len(y1)) )
    verbalise("C", "Initial number of features: %d" % X1.shape[1] )
    verbalise("C", "Final number of features: %d" % ffn )
    estimator = svm.SVR(kernel="linear")
    selector = RFE(estimator, ffn, step=1)
    selector = selector.fit(X1, y1)
    newX1 = selector.transform(X1)
    newX2 = selector.transform(X2)
    return newX1, newX2, selector.support_, selector.ranking_

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

def relativise_vst(X, logpowers=False):
    """
    make each gene expression value relative to mean expression value of all genes. If
    logpowers is False, then the values returned are the power to which the geometric
    mean of the samples **untransformed** values must be raised in order to get the
    untransformed value of the gene. If logpowers is True, then the values are the
    natural log of this power value (ie, positive values are higher than mean, negative
    values are lower than mean).
    """
    powers = X / np.mean(X,0)
    if logpowers:
        logpowers = np.log(powers)
        return logpowers
    else:
        return powers

def multitest(df1, df2, degs=None, name1='X1', name2='X2',
                ffn=100, normtogether=False, normind=False, cv=None):
    "normtogether normalises training and testing data together. "

    verbalise("M", "Predicting %s phase with %s using linear kernel SVM" % (name2, name1))
    if degs:
        dfd = reducedf(df1, degs)
    else:
        dfd = df1
    X1, y1, X2, y2, dfindex = create_xy(dfd, df2, normalise=normtogether)
    if normind:
        X1 = relativise_vst(X1,logpowers=True)
        X2 = relativise_vst(X2,logpowers=True)
    weights, df2_score = predictor_svm(X1, y1, X2, y2)

    #verbalise("M", "Using RFE %s" % (name1))
    #X1s, X2s, rfemask, refweights = get_features_rfe(X1, y1, X2, ffn=ffn)
    #df2_score = predictor_svm(X1s, y1, X2s, y2)

    #verbalise("M", "Using RFECV on %s" % (name1))
    #X1c, X2c, rfecvmask, rfecv_gscores = get_features_rfecv(X1s, y1, X2s, cv=cv)
    #df2_score = predictor_svm(X1c, y1, X2c, y2)

    verbalise("M", "Using Random Forest on %s" % (name1))
    X1, y1, X2, y2, dfindex = create_xy(df1, df2, normalise=normtogether)
    df2_score = random_forest_engine(X1, y1, X2, y2)

    print "-" * 75
    #rfe_features = [ pair[0] for pair in zip(dfindex, rfemask) if pair[1] ]
    #rfecv_features = [ pair[0] for pair in zip(dfindex, rfecvmask) if pair[1] ]
    rfe_features = None
    rfecv_features = None
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

def random_mach(genelist, dftrain, dftest, genenum=5, reps=1000):
    all_scores = []
    for i in range(reps):
        random.shuffle(genelist)
        random_set = genelist[:genenum]
        X1, y1, X2, y2, dfindex = create_xy(dftrain.loc[random_set], dftest)
        X1 = preprocessing.scale(X1)
        X2 = preprocessing.scale(X2)
        weight, score = predictor_svm(X1, y1, X2, y2)
        all_scores.append(score)
    return all_scores

def compare_two(pp, dfpair, training, testing, genelist1, genelist2, name1, name2, genenum=5, reps=1000):
        # run analyses:
        print "analysing gene list 1..."
        genesize1 = len(genelist1)
        all_scores1 = random_mach(genelist1, *dfpair, genenum=genenum, reps=reps)

        print "analysing gene list 2..."
        genesize2 = len(genelist2)
        all_scores2 = random_mach(genelist2, *dfpair, genenum=genenum, reps=reps)

        # create normalised histograms for JSD test:
        hist1 = np.histogram(all_scores1, bins=25, range=(0.,1.))
        hist2 = np.histogram(all_scores2, bins=25, range=(0.,1.))
        norm1 = [ 1.0*val/sum(hist1[0]) for val in hist1[0]]
        norm2 = [ 1.0*val/sum(hist2[0]) for val in hist2[0]]
        jsd = JSD(norm1, norm2)

        # plot histograms:
        plt.hist(all_scores1, bins=25, range=(0.,1.), alpha=0.4, color='red')
        plt.hist(all_scores2, bins=25, range=(0.,1.), alpha=0.4, color='blue')
        title = "Training set: %s | Testing set: %s" % (training, testing)
        plt.title(
            "%s\n%s (red): %d genes | %s (blue): %d genes\nJensen-Shannon divergence: %.3f || KS test: p%s0.05" % (title,
                    name1, genesize1, name2, genesize2,
                    jsd,
                    "<" if ks_2samp(all_scores1, all_scores2)[1]<=0.05 else ">") )
        med1 = np.median(all_scores1)
        med2 = np.median(all_scores2)

        plt.annotate("",
                xy=(med1,0),
                xytext=(med1,max(hist1[0])/10.),
                xycoords='data',
                textcoords='data',
                ha='center',
                arrowprops={'linewidth':1,'arrowstyle':"simple", 'color':'red'})
        plt.annotate("",
                xy=(med2+0.02,0),
                xytext=(med2+0.02,max(hist1[0])/10.),
                xycoords='data',
                textcoords='data',
                ha='center',
                arrowprops={'linewidth':1,'arrowstyle':"simple", 'color':'blue'})

        plt.tight_layout()

        plt.savefig(pp, format='pdf')
        plt.close()

def reducedf(df, genelist):
    #rint df.info()
    #print 'genelist = %s' % genelist
    #print df.loc[genelist].info()
    return df.loc[genelist]

def deg_df_list(deseq_out, orthodf=None, convert=False):
    df = pd.read_csv(deseq_out, header=0, index_col=0, sep=' ')
    try:
        df = df[ df['padj'] <= 0.05 ]
    except KeyError:
        print deseq_out
        print df.columns
        print df
        exit()
    if convert:
        df = df.join(orthodf).dropna()
        degs = df['cbir_loc'].values
    else:
        degs = df.index.values

    print degs[:5]
    return degs

###########  ANALYSES  ################
if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    config_d = load_config(args.configfile[0])

    print "Creating ortholog db..."
    orthodf = create_ortho_df(config_d['convert_orthos'])

    print "Extracting lists of DEGs..."
    # create lists of all DEGs from each experiment:
    meth_degs = deg_df_list(config_d['meth_degs_f'])
    brsw_degs = deg_df_list(config_d['brsw_degs_f'])
    acro_degs = deg_df_list(config_d['acro_degs_f'], orthodf, convert=True)
    bsse_degs = deg_df_list(config_d['bsse_degs_f'])
    sinv_degs = deg_df_list(config_d['sinv_rnr_degs'], orthodf, convert=True)
    cbir_degs = deg_df_list(config_d['cbir_degs_f'])
    dqua_degs = deg_df_list(config_d['dqua_degs_f'], orthodf, convert=True)

    # create dictionary for extracting ncbi names of Cbir genes:
    ncbi = ncbi_dic(file=config_d['ncbi_names'])

    """
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
                for file in os.listdir(config_d['bsse_rd'])
                    if re.search("B[456].*(_FL|_SP).+htseq.gene",file) ]
    sinv_r = [ os.path.join(config_d['sinv_rd'],file)
                for file in os.listdir(config_d['sinv_rd']) if re.search("htseq.gene",file) ]
    """

    #verbalise( "Number of files found for each study:")
    #verbalise("G", "Meth: %d\nBrSw: %d\nBSse: %d\nPolK: %d\nPolU: %d\nAcro: %d" % (len(meth_r), len(brsw_r), len(bsse_r), len(polk_r), len(polu_r), len(acro_r)))

    print "\nAssembling dataframes...\n"
    """
    meth_df = assemble_df(meth_r, lognorm=True, sfnorm=False)
    brsw_df = assemble_df(brsw_r, lognorm=True, sfnorm=False)
    bsse_df = assemble_df(bsse_r, lognorm=True, sfnorm=False)
    polk_df = assemble_df(polk_r, lognorm=True, sfnorm=False)
    polu_df = assemble_df(polu_r, lognorm=True, sfnorm=False)
    acro_df = assemble_df(acro_r, lognorm=True, sfnorm=False)
    cbir_df = meth_df.join(brsw_df).join(bsse_df)
    sinv_df = assemble_df(sinv_r, lognorm=True, sfnorm=False)
    """

    # create dataframes of vst:
    meth_vst = pd.read_csv(config_d['meth_vst'], sep=" ", header=0)
    brsw_vst = pd.read_csv(config_d['brsw_vst'], sep=" ", header=0)
    cbir_vst = pd.read_csv(config_d['cbir_vst'], sep=" ", header=0)
    polu_vst = pd.read_csv(config_d['polu_vst'], sep=" ", header=0)
    polk_vst = pd.read_csv(config_d['polk_vst'], sep=" ", header=0)
    acro_vst = pd.read_csv(config_d['acro_vst'], sep=" ", header=0)
    sinv_vst = pd.read_csv(config_d['sinv_vst'], sep=" ", header=0)
    dqua_vst = pd.read_csv(config_d['dqua_vst'], sep=" ", header=0)

    # convert to orthologs:
    polu_ortho = index_orthos(polu_vst, orthodf)
    polk_ortho = index_orthos(polk_vst, orthodf)
    acro_ortho = index_orthos(acro_vst, orthodf)
    sinv_ortho = index_orthos(sinv_vst, orthodf)
    dqua_ortho = index_orthos(dqua_vst, orthodf)


    all_degs = set( cbir_degs ) | set( acro_degs ) | set( sinv_degs ) | set( dqua_degs )

    # combine all ant orthologs into singe df:
    ants_ortho = cbir_vst.join(acro_ortho).join(sinv_ortho).join(dqua_ortho).dropna()
    verbalise("Y", "Number of genes in common = %-6d total size = %d" % (len(set(ants_ortho.index.values)), len(ants_ortho.index.values)))

    # extract only orthologs that are 1:1:1:1 (some "orthologs" have two genes in
    # some species! (708 to be precise, across all four species after nans have been
    # removed)
    uniq_np, np_counts = np.unique(ants_ortho.index.values, return_counts=True)
    uniq_loci = []
    for i, idx in enumerate(uniq_np):
        if np_counts[i] == 1:
            uniq_loci.append(idx)

    ants_uniq = ants_ortho.loc[uniq_loci,:]

    verbalise("C", "Combined df (4 ants) has %d uniquely orthologous genes and %d samples" % (ants_uniq.shape[0], ants_uniq.shape[1] ))

    ants_by_deg = ants_uniq.loc[list(all_degs)]
    print ants_by_deg.shape
    print len(all_degs), "DEGs with C.biroi ortholog present across all pairwise comparisons"
    print len(cbir_degs), "DEGs for C.biroi"
    print len(sinv_degs), "DEGs for S.invicta"
    print len(set(cbir_degs) | set(sinv_degs)) , "DEGs shared between Cbir and Sinv"
    y = to_binary(ants_uniq.columns)
    X = ants_uniq.T.values
    names = ants_uniq.index.values
    exit()

    # recursive feature elimination with cross validation of combined set:
    """
    verbalise("C", "Sample size of set: %d" % (len(y)) )
    verbalise("C", "Initial number of features: %d" % X.shape[1] )
    cv = 10
    estimator = svm.SVR(kernel="linear")
    selector = RFECV(estimator, cv=cv, step=10, scoring=None)
    selector = selector.fit(X, y)


    reducedX = selector.transform(X)
    print selector.grid_scores_
    verbalise("C",
        "Final number of features: %d (average of %d CV scores = %.3f)" %
            (selector.n_features_, cv, selector.grid_scores_[0] ))

    genescore = zip(selector.support_, selector.ranking_)

    verbalise("C", "%d features:\n%s\nRANKING:\n%s" % (selector.n_features_,
                    [ i for i in range(len(selector.support_)) if selector.support_[i]],
                    len(selector.ranking_)))

    selected_features = [ (show_name(names[i], ncbi), score[1]) for i,score in enumerate(genescore) if selector.support_[i]]
    sorted_features = sorted(selected_features, key=lambda x: x[1], reverse=True)


    verbalise("M", '\n'.join([ "%d %s" % (i[1], i[0]) for i in sorted_features]))


    # plot results:
    plt.figure()
    plt.xlabel("Number of features selected")
    plt.ylabel("Cross validation score (nb of correct classifications)")
    plt.plot(range(1, len(selector.grid_scores_) + 1), selector.grid_scores_)
    plt.show()
    """

    # use to find good number of trees:
    """
    error_rate = []
    for n in np.arange(2,500,2):
        rf = RandomForestClassifier(max_features="log2", n_estimators=n,
                                    oob_score=True,random_state=125)
        #print rf
        rf.fit(X,y)
        error_rate.append( 1 - rf.oob_score_ )

    plt.plot(error_rate)
    plt.xlabel("n_estimators")
    plt.ylabel("OOB error rate")
    plt.show()
    """

    # iteratively drop genes:
    """
    dropem = []
    for i in range(200):
        sh_before = ants_uniq.shape
        df = ants_uniq.drop(dropem)
        sh_after = df.shape
        verbalise("C", sh_before, "minus %d genes =" % len(dropem), sh_after)
        verbalise("B", dropem[-20:])

        # create matrices:
        X = df.T.values
        y = to_binary(df.columns)
        names = df.index.values

        # calculate random forest:
        rf = RandomForestClassifier(max_features='log2', n_estimators=120,
                                    oob_score=True, random_state=125)
        rf.fit(X,y)
        verbalise("Y", "score:", rf.oob_score_)
        # select worst performers to add to dropem list:
        best_to_worst = sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), names),
                                reverse=True)
        verbalise("C", len(best_to_worst))
        dropem += [ p[1] for p in best_to_worst[ -(int(len(best_to_worst)/20.0)): ] ]

        print "\n".join([ "%.4f %s" % (p[0], show_name(p[1],ncbi)) for p in best_to_worst[:5] ])
        print "..."
        print "\n".join([ "%.4f %s" % (p[0], show_name(p[1],ncbi)) for p in best_to_worst[-5:] ])
    """

    #print "Features sorted by their score:"
    #print "\n".join([ "%.4f %s" % (p[0], show_name(p[1],ncbi)) for p in sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), names), reverse=True) if p[0] > 0 ])

    #cum_import=[ f + sum(rf.feature_importances_[:i]) for i,f in enumerate(rf.feature_importances_) if f > 0]
    #print cum_import[:30]

    #plt.plot(cum_import)
    #plt.show()


    # prepare Kegg objects:
    #kegg_tree = kegg.KeggTree(config_d['kegg_f'],config_d['cbir_ko'])

    ###########################  TESTING SITE  #############################

    genelist43 = list(set(["LOC105281428_Q","LOC105278524","LOC105280170","LOC105282036","LOC105275236",
        "LOC105277456","LOC105283266","LOC105284141","LOC105285273","LOC105281340",
        "LOC105286142","LOC105284539","LOC105279365","LOC105279366","LOC105275821",
        "LOC105283812","LOC105287615","LOC105277435","LOC105287013","LOC105282105",
        "LOC105285370","LOC105288082","LOC105287763","LOC105281770","LOC105285920",
        "LOC105281330","LOC105279396","LOC105279397","LOC105285921","LOC105284974",
        "LOC105275269","LOC105285466","LOC105280465","LOC105283721","LOC105275039",
        "LOC105277654","LOC105275039","LOC105277654","LOC105285466","LOC105285597",
        "LOC105288201","LOC105281757","LOC105277039","LOC105284284","LOC105287193",
        "LOC105280785"]))
    #genelist = list(cbir_df.index[:])
    #genelist = list(set(brsw_degs) - (set(meth_degs) & set(brsw_degs)))
    #genelist = list((set(brsw_degs) | set(meth_degs)) - (set(meth_degs) & set(brsw_degs)))
    #genelist = list(set(brsw_degs))
    allcbir = list(set(['LOC105284049','LOC105283429','LOC105279511','LOC105286142',
                'LOC105281037','LOC105279365','LOC105281030','LOC105279366',
                'LOC105287949','LOC105274630','LOC105279939','LOC105283509',
                'LOC105277435','LOC105279564','LOC105281939','LOC105278524',
                'LOC105281330','LOC105276751','LOC105284284','LOC105280537',
                'LOC105276512','LOC105280125','LOC105282643','LOC105275039',
                'LOC105281428_Q','LOC105284181','LOC105276905','LOC105285677',
                'LOC105275269','LOC105275345','LOC105277203','LOC105279847',
                'LOC105274740','LOC105288139','LOC105276234','LOC105285273',
                'LOC105275986','LOC105275513','LOC105287763','LOC105275517',
                'LOC105277039','LOC105279090','LOC105276132','LOC105282911',
                'LOC105276523','LOC105274991','LOC105274638','LOC105279396',
                'LOC105279397','LOC105281340','LOC105280529','LOC105283864',
                'LOC105280581','LOC105279851','LOC105282816','LOC105277793',
                'LOC105286776','LOC105287714','LOC105278588','LOC105278710',
                'LOC105279540','LOC105285144','LOC105285890','LOC105275876',
                'LOC105287013','LOC105276731','LOC105277325','LOC105276925',
                'LOC105275113','LOC105282267','LOC105276628','LOC105276143',
                'LOC105281685','LOC105284141','LOC105279907','LOC105279565',
                'LOC105276495','LOC105288116','LOC105275889','LOC105276092',
                'LOC105286769','LOC105280308','LOC105277547','LOC105287170',
                'LOC105280465','LOC105288064','LOC105279749','LOC105283008',
                'LOC105283644','LOC105280074','LOC105274595','LOC105281444',
                'LOC105280170','LOC105280641','LOC105279871','LOC105280489',
                'LOC105286755','LOC105287193','LOC105274940','LOC105287758',
                'LOC105287615','LOC105281011','LOC105279671','LOC105277260',
                'LOC105278567','LOC105285920','LOC105285376','LOC105286059',
                'LOC105283266','LOC105278238','LOC105277778','LOC105276419',
                'LOC105287790','LOC105278869','LOC105279696','LOC105279164',
                'LOC105277324','LOC105277441','LOC105277876','LOC105277875',
                'LOC105281770','LOC105279766','LOC105279197','LOC105281776',
                'LOC105279590','LOC105281428_W','LOC105284539','LOC105285380',
                'LOC105275762','LOC105279927','LOC105285823','LOC105277349',
                'LOC105287207','LOC105285644','LOC105284533','LOC105283721',
                'LOC105278890','LOC105282036','LOC105281817','LOC105275236',
                'LOC105277456','LOC105287636','LOC105287637','LOC105275428',
                'LOC105285466','LOC105277457','LOC105277458','LOC105287923',
                'LOC105281012','LOC105277690','LOC105281396','LOC105277844',
                'LOC105276346','LOC105282946','LOC105277350','LOC105278469',
                'LOC105274965','LOC105282105','LOC105284285','LOC105283205',
                'LOC105282438','LOC105276066','LOC105285632','LOC105278211',
                'LOC105281708','LOC105281630','LOC105278498','LOC105283374',
                'LOC105276174','LOC105274728','LOC105286727','LOC105276211',
                'LOC105286175','LOC105286889','LOC105280338','LOC105274480',
                'LOC105274936','LOC105280544','LOC105286952','LOC105279728',
                'LOC105285597','LOC105279708','LOC105277621','LOC105286129',
                'LOC105281711','LOC105283337','LOC105280731','LOC105279271',
                        ]))

    ultimate = list(set(["LOC105276751","LOC105279871","LOC105280489","LOC105279696",
                        "LOC105287193","LOC105287636","LOC105281037","LOC105279365",
                        "LOC105277793","LOC105274940","LOC105282816","LOC105287615",
                        "LOC105279564","LOC105282946","LOC105278524","LOC105275876",
                        "LOC105277778","LOC105276512","LOC105282438","LOC105283266",
                        "LOC105276143","LOC105285823","LOC105275269","LOC105276495",
                        "LOC105276234","LOC105279366","LOC105280544","LOC105274480",
                        "LOC105277039","LOC105279197","LOC105281776","LOC105285597",
                        "LOC105274728","LOC105281630","LOC105281711","LOC105287207",
                        "LOC105284533","LOC105280170","LOC105279271",
                        ]))

    combos = [((brsw_vst, meth_vst),
                "brsw PE vst","meth PE vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((cbir_vst,sinv_ortho),
                "all Cbir vst","Sinv vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((dqua_ortho,sinv_ortho),
                "Dqua","Sinv vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((cbir_vst,acro_ortho),
                "all Cbir vst","Acro vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((cbir_vst, dqua_ortho),
                "all Cbir vst","Dqua vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((sinv_ortho, cbir_vst),
                "Sinv vst","all Cbir vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((dqua_ortho, cbir_vst),
                "Dqua vst","all Cbir vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((sinv_ortho, acro_ortho),
                "Sinv vst","acro vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((sinv_ortho, dqua_ortho),
                "Sinv vst","Dqua vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((dqua_ortho, acro_ortho),
                "Dqua vst","Acro vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((acro_ortho, cbir_vst),
                "Acro vst","Cbir vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((acro_ortho, sinv_ortho),
                "Acro vst","Sinv vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),

                ((acro_ortho, dqua_ortho),
                "Acro vst","Dqua vst",
                list(cbir_vst.index[:]),
                list(set(allcbir)),
                "Random","common to all cbir"),


            ]

    # set parameters:
    ####################################################
    #dfpair = ( brsw_df, meth_df )
    #training = "brsw PE scaled"
    #testing = "meth PE scaled"
    #
    #genelist1 = list(cbir_df.index[:])
    #genelist2 = list(set(brsw_degs))
    #
    #name1 = "Random"
    #name2 = "common DEGs(2 exp, 2 anal)"
    ####################################################

    if args.randomsets:
        pp=PdfPages(logfile[:-3] + 'pdf')
        for params in combos:
            verbalise("B", params[1:3],params[5:])
            #pp=PdfPages('Scaled.Mach_learn.%s_%s.%s_%s.pdf' % (testing, training, name1, name2))
            compare_two(pp, *params, genenum=args.genenum, reps=args.replicates)
        pp.close()

    if args.multitest:
        for params in combos:
            verbalise("B", params[1:3],params[5:])
            multitest(params[0][0], params[0][1], degs=params[4],
                        name1=params[1] , name2=params[2], ffn=args.ffn,
                        normtogether=False, normind=True, cv=args.cv)


    ###########################  WINNING ANALYSES  #########################



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

    """
    # Testing differences with Jenson Shannon Divergence
    all_scores1 = []
    for i in range(1000):
        random.shuffle(genelist)
        random_set = genelist[:5]
        X1, y1, X2, y2, dfindex = create_xy(cbir_df.loc[random_set], bsse_df)
        weight, score = predictor_svm(X1, y1, X2, y2)
        all_scores1.append(score)

    all_scores2 = []
    for i in range(1000):
        random.shuffle(genelist)
        random_set = genelist[:45]
        X1, y1, X2, y2, dfindex = create_xy(cbir_df.loc[random_set], bsse_df)
        weight, score = predictor_svm(X1, y1, X2, y2)
        all_scores2.append(score)

    hist1 = np.histogram(all_scores1, bins=24, range=(0.,1.))
    hist2 = np.histogram(all_scores2, bins=24, range=(0.,1.))
    norm1 = [ 1.0*val/sum(hist1[0]) for val in hist1[0]]
    norm2 = [ 1.0*val/sum(hist2[0]) for val in hist2[0]]

    print hist1[1]
    print hist2[1]
    print norm1
    print sum(norm1)
    print sum(norm2)

    verbalise("Y", "Jenson-Shannon Divergence is %.3f" % JSD(norm1, norm2))
    plt.hist(all_scores1, bins=24, range=(0.,1.), alpha=0.4, color='red')
    plt.hist(all_scores2, bins=24, range=(0.,1.), alpha=0.4, color='blue')
    plt.show()


    """


    ########################################################################


