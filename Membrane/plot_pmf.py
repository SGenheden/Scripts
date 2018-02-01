# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to plot PMFs of small molecules pulled through membranes

Very shaky code!
"""

import sys
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pylab as plt
import numpy as np
import os
import argparse

import wham
from sgenlib import colors
from sgenlib import parsing

# The column index into the density arrays
ZCOL=0
WDENS=1
PDENS=2
GDENS=3
TDENS=4

def finalize_plot(fig,min_pmf,max_pmf,filename,ncol=1) :
  #fig.gca().plot([1.0,1.0],[min_pmf-5,max_pmf+5],'k--')
  #fig.gca().plot([1.6,1.6],[min_pmf-5,max_pmf+5],'k--')
  #fig.gca().plot([2.2,2.2],[min_pmf-5,max_pmf+5],'k--')
  fig.gca().set_ylim([min_pmf-5,max_pmf+5])
  #fig.gca().text(0.0,min_pmf-4,"tail",size=8,style='italic')
  #fig.gca().text(1.1,min_pmf-4,"glyc.",size=8,style='italic')
  #fig.gca().text(1.7,min_pmf-4,"head.",size=8,style='italic')
  #fig.gca().text(2.5,min_pmf-4,"water",size=8,style='italic')
  #fig.gca().legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=ncol, mode="expand", borderaxespad=0.)
  fig.gca().legend(loc='best', fancybox=True, fontsize=8,labelspacing=0.20)#framealpha=0.5,
  #fig.gca().set_xlabel("Distance from membrane center [nm]",size=8)
  fig.gca().set_xlabel("Distance from phosphate peak [nm]",size=8)
  fig.gca().set_ylabel("PMF [kJ/mol]",size=8)
  #fig.gca().set_xlim([-32.0,0.0])
  for tick in fig.gca().xaxis.get_major_ticks() :
    tick.label.set_fontsize(8)
  for tick in fig.gca().yaxis.get_major_ticks() :
    tick.label.set_fontsize(8)
#plt.show()
  fig.savefig(filename,format="png",dpi=300)

def _calc_zorg(densities, index):
    di = max(index, len(densities)-1)
    mi = np.argmax(densities[di][:, PDENS])
    return densities[di][mi, ZCOL]

def _plot_areas(axis, densities, centre, ylim):

    def _calc_intercept(density, i, j):
        deltad = density[:,j]/density[:,j].max()-density[:,i]/density[:,i].max()
        return 0.5*(density[np.where(deltad>0)[0],0][-1] +  \
                    density[np.where(deltad<0)[0],0][0])

    y = np.asarray([ylim[0]-5.0,ylim[1]+5.0]*3)
    for i in range(WDENS, GDENS+1):
        x = 0
        for d in densities :
            dhalf = d[d.shape[0]/2:,:]
            xi = _calc_intercept(dhalf, i, i+1)
            if centre == "phosphate" : xi = _calc_zorg(d, i) - xi
            x += xi
        x /= float(len(densities))
        axis.plot([x,x], ylim, 'k--')
    #print x, y




if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Plotting PMFs from umbrella sampling simulations")
    parser.add_argument('-f','--files',nargs="+",help="the output from the calc_pmf.py script.",default=[])
    parser.add_argument('-l','--labels',nargs="+",help="the labels of the PMFs",default=[])
    parser.add_argument('-d','--densities', nargs="+",help="simulation densities",default=[])
    parser.add_argument('-s','--stride',type=float,nargs="+",help="the stride of plot",default=[])
    parser.add_argument('-o','--out',help="the output",default="pmf.png")
    parser.add_argument('-c','--centre',choices=["middle","phosphate"],help="where to centre the PMF",default="middle")
    parser.add_argument('--apl',type=float,help="area per lipid",default=0.67)
    parser.add_argument('--temp',type=float,help="the temperature used in the simulation",default=300.0)
    args = parser.parse_args()

    if len(args.files) == 0 :
        print "No files specified so nothing to do."
        quit()

    boltzmann = 0.001982923700
    RT = args.temp*wham.KB[wham.KJMOL]
    const = args.apl*np.power(10,-24.0)/(786.14*1.66*np.power(10,-27.0))


    densities = []
    if args.densities is not None :
        densities = [parsing.parse2ndarray(d) for d in args.densities]

    f = plt.figure(figsize=(3.33,2.5))
    #a = fig1.add_axes((0.15,0.2,0.8,0.75))
    a = f.add_axes((0.18,0.2,0.75,0.75))

    max_pmf = -1000
    min_pmf = 1000
    pmfs = []
    for i,(filename, label) in enumerate(zip(args.files,args.labels)) :
        pmfs.append(wham.UmbrellaPmf())
        pmfs[-1].read(filename)
        pmfs[-1].z *= 0.1
        transfer_dg = pmfs[-1].transfer_dg()
        wat_barr = pmfs[-1].waterlipid_barrier()
        pen_barr = pmfs[-1].penetration_barrier()
        std_dg = pmfs[-1].standard_dg()

        partstr = ""
        if densities :
            if i < len(densities) :
                part = pmfs[-1].partition(densities[i][:,:WDENS+1])
            else :
                part = pmfs[-1].partition(densities[-1][:,:WDENS+1])
            part_std = part[1]*const
            part     = np.log10(part[0]*const)
            partstr = "| %.3f +- %.3f"%(part,np.abs(part_std/(part*np.log(10))))

        print "%20s %8.3f +- %8.3f | %8.3f +- %8.3f | %8.3f +- %8.3f | %8.3f +- %8.3f | %8.3f at %8.3f %s"%(label,transfer_dg[0],transfer_dg[1],wat_barr[0],wat_barr[1],pen_barr[0],pen_barr[1],std_dg[0],std_dg[1],pmfs[-1].av.min(),pmfs[-1].z[pmfs[-1].av.argmin()],partstr)
        max_pmf = max(max_pmf,pmfs[-1].av.max())
        min_pmf = min(min_pmf,pmfs[-1].av.min())

    for i,(pmf,label) in enumerate(zip(pmfs,args.labels)) :

        if args.centre == "phosphate" and densities:
            pmf.z = _calc_zorg(densities, i) - pmf.z
        kwargs = {}
        if i < len(args.stride) : kwargs["stride"]=args.stride[i]
        pmf.plot(fig=f,label=label,color=colors.color(i),**kwargs)

    #_plot_areas(a, densities, args.centre, (min_pmf-5.0, max_pmf+5.0))

    f.savefig(args.out,format="png",dpi=300)
  #finalize_plot(fig1,min_pmf,max_pmf,args.out,ncol=3)
