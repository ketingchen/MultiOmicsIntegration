import pandas as pd
import numpy as np
from statistics import mean, sqrt
import math
import copy as cp
from scipy import stats
from statsmodels.stats.multitest import multipletests
from scipy.stats import pearsonr

def a_ij(x):
    return(sum(x)-1)

def getWGS(cluster, beta=16, S=500):
    corr0 = abs(cluster.corr(method="pearson"))**beta
    wgs0 = corr0.apply(a_ij, axis=0)
    mu = mean(wgs0)
    WGS = pd.DataFrame()
    for s in range(S):
        idx = np.random.randint(len(list(range(cluster.shape[0]))), size=len(list(range(cluster.shape[0]))))
        idx = idx.tolist()
        cluster_new = cluster.iloc[idx,]
        corr_new = abs(cluster_new.corr(method="pearson"))**beta
        wgs = corr_new.apply(a_ij, axis=0)
        WGS[s] = wgs
        fold = math.floor(S/10)
        if s % fold == 0:
            print(math.floor(s / fold)*10, "%", end=" ->", sep="")
        if s == S-1:
            print("100%")
    return WGS, mu

def getWPlus(WGS):
    ### Iterate by run s out of total runs of S
    mu = WGS[1]
    WGS = WGS[0]
    W_plus = [0] * WGS.shape[0]
    S = [0]*WGS.shape[0]
    X_k = pd.DataFrame(0, index=WGS.index, columns=WGS.columns)
    for s in range(WGS.shape[1]):
        X_k.iloc[:,s] = WGS.iloc[:,s] - mu ### calculate X_k parameter for all genes
    for n in range(WGS.shape[0]):
        xk = [x for x in X_k.iloc[n,:] if x!=0]
        xk2 = [abs(x) for x in xk]
        s = len(xk)
        W = [i[0]+1 for i in sorted(enumerate(xk2), key=lambda x:x[1])]
        idx_aboveZero = [i for i in range(len(xk)) if xk[i]>0]
        W_plus[n] = sum([W[i] for i in idx_aboveZero])
        S[n] = s

    keys = WGS.index.tolist()
    W_plus = dict(zip(keys, W_plus))
    N1 = dict(zip(keys, S))
    return W_plus, N1

def Wtest(WPlus):
    S = WPlus[1]
    WPlus = WPlus[0]
    keys = [*WPlus]
    p_values = [1] * len(WPlus)
    i = 0
    for k in keys:
        wp = WPlus.get(k)
        s = S.get(k)
        E_H0_Wplus = s * (s + 1) / 4
        V_H0_Wplus = s * (s + 1) * (2 * s + 1) / 24
        test_statistic = (wp - E_H0_Wplus) / sqrt(V_H0_Wplus)
        p_values[i] = 1 - stats.norm(0, 1).cdf(test_statistic)
        i += 1
    p_adj = multipletests(p_values)[1].tolist()
    result = pd.DataFrame({'Gene':keys, 'p_values':p_values, 'p_adj':p_adj})
    return result

def geneModuleMembership(genes, module):
    corr = [0] * genes.shape[1]
    pval = [0] * genes.shape[1]
    for i in range(genes.shape[1]):
        result = pearsonr(genes.iloc[:,i], module)
        corr[i] = result[0]
        pval[i] = result[1]
    return pd.DataFrame({'Gene':genes.columns, 'corr':corr, 'pval':pval})

def geneTraitSignificance(genes, trait):
    corr = [0] * genes.shape[1]
    pval = [0] * genes.shape[1]
    for i in range(genes.shape[1]):
        result = pearsonr(genes.iloc[:,i], trait)
        corr[i] = result[0]
        pval[i] = result[1]
    return pd.DataFrame({'Gene':genes.columns, 'corr':corr, 'pval':pval})

def main(cluster, beta=16, S=1000):
    WGS = getWGS(cluster, beta, S)
    WPlus = getWPlus(WGS)
    p_values = Wtest(WPlus)
    p_values = p_values.sort_values('p_adj')
    return p_values

def findHub(p_values, cluster, traits, module, corr_cut=0.5, whichTrait="both"):
    cls_sig = cluster[[p_values.iloc[i,0] for i in range(p_values.shape[0]) if p_values.iloc[i,2]==0]]
    gm_corr = geneModuleMembership(cls_sig, module)
    gt_corr_cutin = geneTraitSignificance(cls_sig, traits['cutin'])
    gt_corr_ecw = geneTraitSignificance(cls_sig, traits['ecw'])
    if whichTrait == "both":
        return gm_corr, gt_corr_cutin, gt_corr_ecw, [gm_corr.iloc[i,0] for i in range(gm_corr.shape[0]) if gm_corr.iloc[i,1]>0.8 and abs(gt_corr_cutin.iloc[i,1])>corr_cut and abs(gt_corr_ecw.iloc[i,1])>corr_cut]
    elif whichTrait == "cutin":
        return gm_corr, gt_corr_cutin, gt_corr_ecw, [gm_corr.iloc[i,0] for i in range(gm_corr.shape[0]) if gm_corr.iloc[i,1]>0.8 and abs(gt_corr_cutin.iloc[i,1])>corr_cut]
    else:
        return gm_corr, gt_corr_cutin, gt_corr_ecw, [gm_corr.iloc[i,0] for i in range(gm_corr.shape[0]) if gm_corr.iloc[i,1]>0.8 and abs(gt_corr_ecw.iloc[i,1])>corr_cut]

def edgelist(cluster, beta=16, mu=0):
    corr0_beta1 = cluster.corr(method="pearson")
    corr0 = corr0_beta1 ** beta
    corr0_abs = abs(corr0)
    wgs0 = corr0_abs.apply(a_ij, axis=0)
    if mu == 0:
        mu = sum(wgs0) / (corr0.shape[0]**2 - corr0.shape[0])
    corr1 = corr0.mask(np.tril(np.ones(corr0.shape)).astype(np.bool))
    corr1 = corr1.stack().reset_index()
    corr1.columns = ['Gene_A', 'Gene_B', 'corr']
    corr1_beta1 = corr0_beta1.mask(np.tril(np.ones(corr0_beta1.shape)).astype(np.bool))
    corr1_beta1 = corr1_beta1.stack().reset_index()
    mask = (abs(corr1['corr']) > mu)
    edl = cp.deepcopy(corr1[mask])
    edl['sign'] = corr1_beta1.iloc[:,2]>0
    return edl