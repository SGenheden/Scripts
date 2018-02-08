# Author: Samuel Genheden samuel.genheden@gmail.com

import sys
import os

import numpy as np
import matplotlib.pylab as plt
import scipy.stats as stats

from sgenlib import colors
from sgenlib import parsing

def _plot_data(figure, datalist, labels, ylabels, xlabels, ncols=3):
    if isinstance(ylabels,str):
        ylabels = [ylabels]*len(datalist)
        xlabels = [xlabels]*len(datalist)
        ilabels = False
    else:
        ilabels = True

    nrows = int(np.ceil(len(datalist)/float(ncols)))
    minv = np.floor(min([d.min() for d in datalist]))
    maxv = np.ceil(max([d.max() for d in datalist]))
    vrange = maxv - minv
    delta = np.ceil(vrange/8.0)
    ticks = np.arange(minv,maxv+1.0,delta)

    for i,(data,label,ylabel,xlabel) in enumerate(zip(datalist,labels,ylabels,xlabels),1):
        a = figure.add_subplot(nrows,ncols,i)
        a.plot(data[:,1],data[:,0],".",color=colors.color(i-1))
        tau,tpval = stats.kendalltau(data[:,1],data[:,0])
        r,rpval = stats.pearsonr(data[:,1],data[:,0])
        print "%s\t%.5f\t%.5f\t%.5f\t%5f"%(label,r,rpval,tau,tpval)
        a.plot([minv,maxv],[minv,maxv],"--k")
        a.set_xlim(minv,maxv)
        a.set_ylim(minv,maxv)
        a.set_xticks(ticks)
        a.set_yticks(ticks)
        if ilabels or a.is_first_col():
            a.set_ylabel(ylabel)
        if ilabels or a.is_last_row():
            a.set_xlabel(xlabel)
        a.text(0.05,0.92,label,transform=a.transAxes)
        a.set_aspect('equal')

    figure.tight_layout()
    print ""

if __name__ == '__main__' :

    data = parsing.parse2ndarray(sys.argv[1])
    figure = plt.figure(1,figsize=(7,7))
    _plot_data(figure, [data], [""],r'$\mathrm{log} D_\mathrm{AA}$',r'$\mathrm{log} D_\mathrm{exp.}$',ncols=1)
    figure.savefig(sys.argv[2], format="png", dpi=300)
