# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse

import numpy as np
import matplotlib.pylab as plt

import sheetslib
from sgenlib import colors

def _make_plot(residue, analytical, observed, error, axis):

    left1 = np.arange(0.4,0.4+2.0*len(residue),2.0)
    left2 = left1 + 0.8

    abar = axis.bar(left1, analytical, width=0.8, color='w')
    obar = axis.bar(left2, observed, yerr=error, ecolor="k", width=0.8, color='w', hatch='+')
    for l, a, o, e in zip(left2, analytical, observed, error):
        if np.abs(a-o) > 1.96*e and o > a :
            axis.text(l+0.2,o+0.05,"*")

    axis.set_xticks(left2)
    axis.set_xticklabels(residue)
    axis.set_ylim([0,1.0])
    for tick in axis.xaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    for tick in axis.yaxis.get_major_ticks() :
        tick.label.set_fontsize(8)

    return abar, obar

def _extract_data(sheet):

    residues = []
    analyticals = []
    observeds = []
    errors = []
    maxres = 0
    for sysoff in [0,7,14]:
        res = sheetslib._extract_residues(sheet,rowstart=2,sysoffset=sysoff)
        maxres = max(maxres,len(res))

    for sysoff in [0,7,14]:
        res = sheetslib._extract_residues(sheet,rowstart=2,sysoffset=sysoff)
        nres = len(res)
        res = [res[i][0:3] if i < nres else "" for i in range(maxres)]
        residues.append(res)
        analyticals.append(np.asarray([sheet.cell(column=sysoff+3,row=i+3).value if i < nres else 0.0
                            for i in range(maxres)]))
        observeds.append(np.asarray([sheet.cell(column=sysoff+4,row=i+3).value if i < nres else 0.0
                            for i in range(maxres)]))
        errors.append(np.asarray([sheet.cell(column=sysoff+5,row=i+3).value if i < nres else 0.0
                            for i in range(maxres)]))
    return residues, analyticals, observeds, errors

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to make bar plots")
    argparser.add_argument('-xls','--xls',help="the filename of the XLSX file")
    args = argparser.parse_args()

    wb = sheetslib.open_book(args.xls)
    residues, analyticals, observeds, errors = _extract_data(wb["Bedroc_var"])

    props = dict(facecolor='white',edgecolor='white', alpha=0.8)
    letters = ["A)","B)","C)","D)"]
    labels = ["NH","Met","Dict"]

    f = plt.figure(1,figsize=(7,6))
    for i,(residue, analytical, observed, error) in enumerate(zip(residues, analyticals, observeds, errors),1):
        a = f.add_subplot(3,1,i)
        abar, obar = _make_plot(residue, analytical, observed, error, a)
        #a.text(0.05,0.04,labels[i-1],transform=a.transAxes,bbox=props, size=10)
        a.text(0.03,0.88,letters[i-1]+" "+labels[i-1],transform=a.transAxes,bbox=props, size=10)
        if i == 2:
            a.legend([abar, obar],["Analytical","Observed"],
                    loc=1, fancybox=True, framealpha=0.5,fontsize=8,labelspacing=0.20)
            a.set_ylabel("BEDROC value")
    f.savefig("bedroc_var.png",dpi=300)
