# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse
import os

import openpyxl as xl
import numpy as np
import scipy.stats as stats

import sheetslib
import quality

def _make_array(data1, data2, error1, error2, residues):

    array = []
    for res in residues :
        if res in data1 and res in data2 :
            array.append([data1[res],error1[res],data2[res],error2[res]])
    return array

def _parse_file(filename):

    with open(filename,"r") as f :
        lines = f.readlines()
    rawdata = [line.strip().split() for line in lines]
    data = {}
    errors = {}
    for i, fdata in enumerate(rawdata,2):
        residue = fdata[0]+fdata[1]+"-%d"%i
        data[residue] = float(fdata[2])
        errors[residue] = float(fdata[3])
    return data, errors

def _bootstrap_data(data, metric_func, nboots=500):
    boots = np.zeros(nboots)
    for i in range(nboots):
        idx = np.random.randint(0,data.shape[0]-1,data.shape[0])
        boots[i] = metric_func(data[idx,:])
    return boots.std()

def _mad(data):
    return np.mean(np.abs(data[:,0]-data[:,2]))

def _mud(data):
    return np.median(np.abs(data[:,0]-data[:,2]))

def _msd(data):
    return np.mean(data[:,0]-data[:,2])

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to compared methods")
    argparser.add_argument('-xls','--xls',help="the filename of the XLSX file")
    argparser.add_argument('-sheet','--sheet',help="the sheet in the XLSX file")
    #argparser.add_argument('-sys','--sys',choices=["bpti","lac","l02"],help="the system")
    argparser.add_argument('-metric','--metric',choices=["aud","mud","msd"],help="the metric", default="aud")
    args = argparser.parse_args()

    sysoff = {"bpti":0,"lac":12,"l02":22}
    ffs = ["tip","elba","gb"]
    freqs = [500, 1000, 2500, 5000, 10000, 25000, 50000, 100000, 250000]
    metric_func = {"aud":_mad, "mud":_mud, "msd":_msd}

    wb = sheetslib.open_book(args.xls)

    all_arr = {ff : {f : [] for f in freqs} for ff in ffs}
    for sys in ["bpti","lac","l02"] :
        res, labels, data, errors = sheetslib.extract_sheet(wb, args.sheet, sysoff[sys])
        for freq in freqs :
            print "%d"%(freq/1000),
            for ff in ffs:
                if sys == "lac" and ff == "gb" :
                    continue
                computed, comp_err = _parse_file("s2_nh_%s_%sav_f%d.out"%(sys,ff,freq))
                arr = _make_array(data[0], computed,  errors[0], comp_err, res)
                if sys != "bpti" :
                    all_arr[ff][freq].extend(arr)
                arr = np.asarray(arr)
                print "\t%.5f\t%.5f"%(metric_func[args.metric](arr),_bootstrap_data(arr, metric_func[args.metric])),
            print ""
        print "---"

    for freq in freqs :
        print "%d"%(freq/1000),
        for ff in ffs:
            #print len(all_arr[ff][freq])
            arr = np.asarray(all_arr[ff][freq])
            print "\t%.5f\t%.5f"%(metric_func[args.metric](arr),_bootstrap_data(arr, metric_func[args.metric])),
        print ""

    print "N="
    for ff in ffs:
        print len(all_arr[ff][freqs[0]])
