# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse
import os

import openpyxl as xl
import numpy as np
import scipy.stats as stats

import sheetslib
import quality

def _extract_stride(filename, offset):

    if filename is None : return []

    secondary = []
    with open(filename, "r") as f:
        for line in f.readlines():
            if line[:3] != "ASG" : continue
            struct = line[24]
            if struct in ["H","G","I","E"] :
                secondary.append("%s%d"%(line[5:8],int(line[11:15])+offset))
    return secondary

def _make_array(data1, data2, error1, error2, residues, secondary):

    array = []
    for res in residues :
        if res in data1 and res in data2 and (not secondary or res.split("-")[0] in secondary) :
            array.append([data1[res],error1[res],data2[res],error2[res]])
    return np.asarray(array)

def _calc_mad(data1,data2,residues,secondary,logit=False):
    data = _make_array(data1,data2,residues,secondary,logit)
    return np.abs(data[:,0]-data[:,1]).mean()

def _calc_median(data1,data2,residues,secondary,logit=False):
    data = _make_array(data1,data2,residues,secondary,logit)
    return np.median(np.abs(data[:,0]-data[:,1]))

def _calc_max(data1,data2,residues,secondary,logit=False):
    data = _make_array(data1,data2,residues,secondary,logit)
    return np.abs(data[:,0]-data[:,1]).max()

def _calc_r(data1,data2,residues,secondary,logit=False):
    data = _make_array(data1,data2,residues,secondary,logit)
    return stats.pearsonr(data[:,0],data[:,1])[0]

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to compared methods")
    argparser.add_argument('-xls','--xls',help="the filename of the XLSX file")
    argparser.add_argument('-sheet','--sheet',help="the sheet in the XLSX file")
    argparser.add_argument('-stride','--stride',help="the stride output")
    argparser.add_argument('-strideoff','--strideoff',type=int,help="the stride output offset",default=0)
    argparser.add_argument('-sysoff','--sysoff',type=int,help="the system column offset",default=0)
    args = argparser.parse_args()

    wb = sheetslib.open_book(args.xls)
    res, labels, data, errors = sheetslib.extract_sheet(wb, args.sheet, args.sysoff)

    secondary = _extract_stride(args.stride,args.strideoff)
    #print secondary

    mads = []
    medians = []
    mses = []
    for i in range(len(labels)):
        for j in range(i+1,len(labels)):
            seldata = _make_array(data[i], data[j],
                                    errors[i], errors[j], res, secondary)
            qc = quality.QualityCollection(seldata[:,0],seldata[:,1],
                                            seldata[:,2],seldata[:,3],
                                            ["mad","absmed","msd"],nboots=5)

            qc.bootstrap()
            mads.append((qc.metriclst[0].biased,qc.metriclst[0].std))
            medians.append((qc.metriclst[1].biased,qc.metriclst[1].std))
            mses.append((qc.metriclst[2].biased,qc.metriclst[2].std))

    print "\t",
    for i in range(len(labels)):
        for j in range(i+1,len(labels)):
            print "%s-%s\t\t"%(labels[i],labels[j]),
    print ""

    print "MAD\t",
    print "\t".join("%.4f\t%.4f"%v for v in mads)
    print "MED\t",
    print "\t".join("%.4f\t%.4f"%v for v in medians)
    print "MSE\t",
    print "\t".join("%.4f\t%.4f"%v for v in mses)
