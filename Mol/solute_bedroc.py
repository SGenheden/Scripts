# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to compute BEDROC metrics for a set of solutes

The -f file should contain the true value in the first
column and the test value in the second column

For each super group, it will compute:
    - BEDROC metric assuming a uniform distribution
    - BEDROC metric on the observed data
    - A p-value on the difference between uniform and observed BEDROC
    - MAD and a p-value comparing to the total MAD
    - MSD

Examples:
solute_bedroc.py -f results.dat -l labels.dat -g Groups/
"""

import argparse
import os

import numpy as np
import croc
import scipy.stats as stats

from sgenlib import parsing

_supergroups = ['alcohol',
                'aldehyde',
                'alkane',
                'alkene',
                'alkyl bromide',
                'alkyl chloride',
                'amine',
                'aromatic compound',
                'aromatic amine',
                'carbonitrile',
                'carboxylic acid',
                'carboxylic acid ester',
                'carboxylic acid amide',
                'ether',
                'halogen derivative',
                'heterocyclic compound',
                'ketone',
                'nitro compound',
                'oxo(het)arene',
                'phenol or hydroxyhetarene']

def _transform_groups(groups):
    """Transform a set of groups to a list of supergroups"""
    groups2 = []
    for group in groups:
        group2 = group.replace("primary","")
        group2 = group2.replace("secondary","")
        group2 = group2.replace("tertiary","")
        if group2.find("ether") > -1:
            group2 = "ether"
        elif group2.find("aliphatic amine") > -1:
            group2 = "amine"
        group2 = group2.strip()
        if group2 in _supergroups:
            groups2.append(group2)
        if group2.find("fluoride") > -1 or group2.find("chloride") > -1 or group2.find("bromide") > -1:
            groups2.append("halogen derivative")
        if group2.find("aromatic amine") > -1:
            groups2.append("aromatic amine")
    return list(set(groups2))

def _bedroc_analytic(npos, ntot):
    npos = npos / float(ntot)
    nneg = 1.0 - npos
    num1 = np.exp(npos)-nneg
    denom1 = np.exp(npos)-1
    denom2 = 1 - np.exp(-nneg)
    return (num1/denom1)-(nneg / denom2)

def _bootstrap_bedroc(data,nboots=500):
    bedrocs = np.zeros(nboots)
    for i in range(nboots):
        while True:
            idx = np.random.randint(0,data.shape[0]-1,data.shape[0])
            if np.sum(data[idx,1]) > 0 and np.sum(data[idx,1]) != data.shape[0]: break
        scoreddata = croc.ScoredData(data[idx,:])
        bedrocs[i] = croc.BEDROC(scoreddata,1.0)['BEDROC']
    return bedrocs

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to compute BEDROC analysis for a groups of solutes")
    argparser.add_argument('-f', '--file', help="the results")
    argparser.add_argument('-l', '--labels', help="the labels")
    argparser.add_argument('-g', '--groupfolder', help="the folder for groups")
    argparser.add_argument('--debug', action="store_true", help="turn on debugging", default=False)
    argparser.add_argument('--plot', action="store_true", help="turn on plotting", default=False)
    args = argparser.parse_args()

    data = parsing.parse2ndarray(args.file)
    with open(args.labels, 'r') as f:
        labels = [line.strip() for line in f.readlines()]

    solutegroups = []
    for label in labels:
        with open(os.path.join(args.groupfolder, label+".groups"), 'r') as f:
            groups = [line.strip() for line in f.readlines()]
        solutegroups.append(groups)

    if args.debug:
        print "Unique basic groups and their counts"
        uniquelist = {}
        for g in solutegroups:
            for l in g :
                if l not in uniquelist : uniquelist[l] = 0
                uniquelist[l] += 1
        for g in uniquelist :
            if _transform_groups([g]) and _transform_groups([g])[0] in _supergroups: print "* ",
            print g, uniquelist[g]
        print "-----\n"        

    errors = data[:,1] - data[:,0]
    abserrors = np.abs(errors)
    nsolutes = errors.shape[0]
    print "Group\tN\tBEDROC uni\tBEDROC obs.\tstd\tsign.\tMSD\tp-value\tMSE"
    for group in _supergroups:
        scoreddata = croc.ScoredData()
        bootdata = np.zeros((nsolutes,2))
        for i, (abserror, solgroup) in enumerate(zip(abserrors, solutegroups)):
            if group in _transform_groups(solgroup) :
                scoreddata.add(abserror, 1)
                bootdata[i,:] = [abserror, 1]
            else:
                scoreddata.add(abserror, 0)
                bootdata[i,:] = [abserror, 0]
        npos = np.sum(bootdata[:,1])
        if npos >= 5 and npos <= nsolutes-5:
            bedroc0 = croc.BEDROC(scoreddata,1.0)['BEDROC']
            bedroc_anal = _bedroc_analytic(npos, nsolutes)
            bedroc_boot = _bootstrap_bedroc(bootdata)
            signstr = " "
            if bedroc_boot.mean() > bedroc_anal :
                t = (bedroc_boot.mean()-bedroc_anal)/bedroc_boot.std()
                pval = stats.t.sf(np.abs(t), npos-1) * 2
                signstr = "%.5f"%pval
            pos_abserr = abserrors[bootdata[:,1]==1]
            pos_err = errors[bootdata[:,1]==1]
            t,pval = stats.ttest_ind(pos_abserr, abserrors, equal_var=False)
            print "%s\t%d\t%.3f\t%.3f\t%.3f\t%s\t%.3f\t%.3f\t%.3f"%(group,
                    npos, bedroc_anal, bedroc_boot.mean(), bedroc_boot.std(),
                    signstr, pos_abserr.mean(), pval, pos_err.mean())
