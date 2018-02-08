# Author: Samuel Genheden samuel.genheden@gmail.com
"""
This is a program to analyse error of order parameters and populate a spreadsheet
"""

import argparse
import os
from collections import namedtuple

import numpy as np
import croc

import sheetslib

OrderParam = namedtuple("OrderParam",["resname","resid","value"])

def _average_group(paramlist, residues, hasexp) :

    avparams = []
    nfilesinv = 1.0 / float(len(paramlist))
    #paramlist = paramlist[::-1]
    offset = 1 if hasexp else 0
    for res in residues:
        if not hasexp :
            av = sum([f[res] for f in paramlist]) * nfilesinv
        for f in paramlist[offset:]:
            if res[:3] in ["HID","HIE","HIP"]:
                resnam = "HIS"
            elif res[:3] in ["CYX","CYS"]:
                resnam = "CYS"
            else:
                resnam = res[:3]
            if hasexp and res in paramlist[0] :
                avparams.append(OrderParam(resnam,res[3:],abs(f[res]-paramlist[0][res])))
            elif not hasexp :
                avparams.append(OrderParam(resnam,res[3:],abs(f[res]-av)))
    return avparams

def _bedroc_analytic(npos,ntot):
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
    return bedrocs.mean(),bedrocs.std()

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Program to perform error analysis of order parameters")
    argparser.add_argument('-xls','--xls',help="the filename of the XLSX file")
    argparser.add_argument('-sheet','--sheet',help="the sheet in the XLSX file")
    argparser.add_argument('-sysoffs','--sysoffs',nargs="+",type=int,help="the system column offset")
    argparser.add_argument('--hasexp',action="store_true",default=False)
    args = argparser.parse_args()

    wb = sheetslib.open_book(args.xls)

    allparams = []
    for sysoff in args.sysoffs:
        res, labels, data, errors = sheetslib.extract_sheet(wb, args.sheet, sysoff)
        avparams = _average_group(data, res, args.hasexp)
        allparams.extend(avparams)

    resnames = sorted(list(set([p.resname for p in allparams])))
    for resname in resnames:
        scoreddata = croc.ScoredData()
        bootdata = np.zeros((len(allparams),2))
        for i,param in enumerate(allparams):
            if param.resname == resname :
                scoreddata.add(param.value,1)
                bootdata[i,:] = [param.value,1]
            else:
                scoreddata.add(param.value,0)
                bootdata[i,:] = [param.value,0]
        if np.sum(bootdata[:,1]) >= 5:
            bedroc_anal = _bedroc_analytic(np.sum(bootdata[:,1]),len(allparams))
            bedroc0 = croc.BEDROC(scoreddata,1.0)['BEDROC']
            bedroc_mean,bedroc_std = _bootstrap_bedroc(bootdata)
            print "%s\t%d\t%.3f\t%.3f\t%.3f"%(resname,np.sum(bootdata[:,1]),
                            bedroc_anal,bedroc0,bedroc_std)
