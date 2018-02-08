# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse

import numpy as np
import matplotlib.pylab as plt

import sheetslib
from sgenlib import colors

def _make_plot(data, errors, axis):

    left = np.asarray([0.8,1.6,2.4,4.0,4.8,5.6,7.2,8.0,8.8])
    #color = [colors.color(0),colors.color(1),colors.color(2)]*3
    color = [(255.0/255.0,255.0/255.0,255.0/255.0),
            (127.0/255.0,127.0/255.0,127.0/255.0),
            (50.0/255.0,50.0/255.0,50.0/255.0)]*3

    lbls = ["AA","CG","GB"]*3
    lbls[5] = ""
    lbls[1] += "\n$Bpti$"
    lbls[4] += "\n$Gal-lac$"
    lbls[7] += "\n$Gal-l02$"

    axis.bar(left,data,width=0.8, yerr=errors, ecolor="k", color=color)
    for l, aa, cg, eaa, ecg in zip(left[1::3], data[::3], data[1::3], errors[::3], errors[1::3]):
        if cg > 0 and np.abs(aa-cg) > 1.96*np.sqrt(eaa**2+ecg**2):
            axis.text(l+0.3,cg+0.05+ecg,"*")
    for l, aa, gb, eaa, egb in zip(left[2::3], data[::3], data[2::3], errors[::3], errors[2::3]):
        if gb > 0 and np.abs(aa-gb) > 1.96*np.sqrt(eaa**2+egb**2):
            axis.text(l+0.3,gb+0.05+egb,"*")
    axis.set_xticks(left+0.4)
    axis.set_xticklabels(lbls)
    axis.set_ylabel("Hydrogen bonds")

    for tick in axis.xaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    for tick in axis.yaxis.get_major_ticks() :
        tick.label.set_fontsize(8)

def _extract_data(sheet):

    def _stde(value):
        if isinstance(value,float) or isinstance(value,int) :
            return value

        stdstr, nstr = value.split("/SQRT")
        std = float(stdstr[1:])
        n = float(nstr[1:-1])
        return std / np.sqrt(n)

    data = []
    errors = []
    for col in [3,5,7]:
        data.append(
            np.asarray([sheet.cell(column=col,row=i).value for i in range(2,11)])
            )
        errors.append(
            np.asarray([_stde(sheet.cell(column=col+1,row=i).value) for i in range(2,11)])
            )

    data.append(np.asarray([data[0][i]-data[1][i]-data[2][i] for i in range(data[0].shape[0])]))
    errors.append(np.asarray([np.sqrt(errors[0][i]**2+errors[1][i]**2+errors[2][i]**2)
                    for i in range(data[0].shape[0])]))

    return data[1:], errors[1:]

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to make bar plots")
    argparser.add_argument('-xls','--xls',help="the filename of the XLSX file")
    argparser.add_argument('-out','--out',help="the output prefix")
    args = argparser.parse_args()

    wb = sheetslib.open_book(args.xls)
    data, errors = _extract_data(wb["Hbond"])

    props = dict(facecolor='white',edgecolor='white', alpha=0.8)
    letters = ["A) BB-BB","B) SC-SC","C) BB-SC"]
    widths = [0,3.25,7]

    f = plt.figure(1,figsize=(7,5))

    for i,(d, e, l) in enumerate(zip(data, errors, letters),1):
        a = f.add_subplot(2,2,i)
        _make_plot(d, e, a)
        #a.text(0.05,0.04,labels[i-1],transform=a.transAxes,bbox=props, size=10)
        a.text(0.03,0.88,letters[i-1],transform=a.transAxes,bbox=props, size=10)
        #if i == 2:
        #    a.legend([abar, obar],["Analytical","Observed"],
        #            loc=1, fancybox=True, framealpha=0.5,fontsize=8,labelspacing=0.20)
        #    a.set_ylabel("BEDROC value")
    f.savefig("hbond.png",dpi=300)
