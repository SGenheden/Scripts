# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse

import numpy as np
import matplotlib.pylab as plt

import sheetslib
from sgenlib import colors

def _make_plot(data, lbls, axis):

    left = np.asarray([0.8,1.6,2.4,4.0,4.8,5.6,7.2,8.0,8.8])
    #color = [colors.color(0),colors.color(1),colors.color(2)]*3
    color = [(255.0/255.0,255.0/255.0,255.0/255.0),
            (127.0/255.0,127.0/255.0,127.0/255.0),
            (50.0/255.0,50.0/255.0,50.0/255.0)]*3

    lbls2 = list(lbls)
    if data[1] != 0 :
        lbls2[1] += "\n$Bpti$"
    lbls2[4] += "\n$Gal-lac$"
    lbls2[7] += "\n$Gal-l02$"

    axis.bar(left,data,width=0.8,color=color)
    axis.set_xticks(left+0.4)
    axis.set_xticklabels(lbls2)

    for tick in a.xaxis.get_major_ticks() :
        tick.label.set_fontsize(8)
    for tick in a.yaxis.get_major_ticks() :
        tick.label.set_fontsize(8)

def _extract_data(sheet, rowoffset):

    colstart = 35 # 44 for secondary structure
    rowstarts = [4,9,14]
    #rowoffset = 2 # 1 for MED and 2 for MSE

    data = []
    labels = []
    for rowstart in rowstarts:
        for coli in range(1,6,2):
            cell = sheet.cell(column=colstart+coli,row=rowstart+rowoffset)
            if cell.value is None :
                data.append(0.0)
                labels.append("")
            else:
                data.append(cell.value)
                labels.append(
                    sheet.cell(column=colstart+coli,row=rowstart-1).value.split("-")[0])

    data = np.asarray(data)
    return data,labels

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to make bar plots")
    argparser.add_argument('-xls','--xls',help="the filename of the XLSX file")
    argparser.add_argument('-sheets','--sheets',nargs="+",help="the sheet in the XLSX file")
    argparser.add_argument('-labels','--labels',nargs="+",help="the labels")
    argparser.add_argument('-out','--out',help="the output prefix")
    argparser.add_argument('-yunit','--yunit',help="the unit of the y axis",default="")
    args = argparser.parse_args()

    wb = sheetslib.open_book(args.xls)

    if len(args.sheets) == 2 :
        ncol = 1
        nrow = 2
    else:
        ncol = 2
        nrow = int(np.ceil(len(args.sheets)/2.0))

    props = dict(facecolor='white',edgecolor='white', alpha=0.8)
    letters = ["A)","B)","C)","D)"]
    widths = [0,3.25,7]

    figs = [plt.figure(i,figsize=(widths[ncol],5))
            for i in range(1,4)]
    ylabels = ["AUD","MUD","MSD"]
    for i,(sheet,label) in enumerate(zip(args.sheets,args.labels),1):
        ws = None
        try :
            ws = wb[sheet]
        except :
            raise Exception("Could not find sheet %s"%sheet)
        for rowoffset,(f,ylabel) in enumerate(zip(figs,ylabels)):
            data, labels = _extract_data(ws, rowoffset)
            a = f.add_subplot(nrow,ncol,i)
            _make_plot(data, labels, a)
            a.text(0.05,0.04,label,transform=a.transAxes,bbox=props, size=10)
            a.text(0.03,0.92,letters[i-1],transform=a.transAxes,bbox=props, size=10)
            a.set_ylabel(ylabel+args.yunit)
            print ylabel+args.yunit
    for f, lbl in zip(figs,ylabels):
        f.tight_layout()
        f.savefig("%s-%s.png"%(args.out,lbl.lower()),format="png",dpi=300)
