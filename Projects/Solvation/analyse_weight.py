# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the weights of all solutes in a list
The solutes are taken from the Minnesota solvation database
"""

import argparse
import re

import numpy as np
import matplotlib.pylab as plt

import dblib

mass = { 'H'  :   1.0079 , 'He' :   4.0026 , 'Li' :   6.941  ,
         'Be' :   9.0122 , 'B'  :  10.811  , 'C'  :  12.0107 ,
         'N'  :  14.0067 , 'O'  :  15.9994 , 'F'  :  18.9984 ,
         'Ne' :  20.1797 , 'Na' :  22.9898 , 'Mg' :  24.3050 ,
         'Al' :  26.9815 , 'Si' :  28.0855 , 'P'  :  30.9738 ,
         'S'  :  32.065  , 'Cl' :  35.453  , 'Ar' :  39.948  ,
         'K'  :  39.0983 , 'Ca' :  40.078  , 'Sc' :  44.9559 ,
         'Ti' :  47.867  , 'V'  :  50.9415 , 'Cr' :  51.9961 ,
         'Mn' :  54.9380 , 'Fe' :  55.845  , 'Co' :  58.9331 ,
         'Ni' :  58.6934 , 'Cu' :  63.546  , 'Zn' :  65.409  ,
         'Ga' :  69.723  , 'Ge' :  72.64   , 'As' :  74.9216 ,
         'Se' :  78.96   , 'Br' :  79.904  , 'Kr' :  83.798  ,
         'Rb' :  85.4678 , 'Sr' :  87.62   , 'Y'  :  88.9059 ,
         'Zr' :  91.224  , 'Nb' :  92.9064 , 'Mo' :  95.94   ,
         'Tc' :  98.     , 'Ru' : 101.07   , 'Rh' : 102.9055 ,
         'Pd' : 106.42   , 'Ag' : 107.8682 , 'Cd' : 112.411  ,
         'In' : 114.818  , 'Sn' : 118.710  , 'Sb' : 121.760  ,
         'Te' : 127.60   , 'I'  : 126.9045 , 'Xe' : 131.293  ,
         'Cs' : 132.9055 , 'Ba' : 137.327  , 'La' : 138.9055 ,
         'Ce' : 140.116  , 'Pr' : 140.9077 , 'Nd' : 144.242  ,
         'Pm' : 145.     , 'Sm' : 150.36   , 'Eu' : 151.964  ,
         'Gd' : 157.25   , 'Tb' : 158.9254 , 'Dy' : 162.500  ,
         'Ho' : 164.9303 , 'Er' : 167.259  , 'Tm' : 168.9342 ,
         'Yb' : 173.04   , 'Lu' : 174.967  , 'Hf' : 178.49   ,
         'Ta' : 180.9479 , 'W'  : 183.84   , 'Re' : 186.207  ,
         'Os' : 190.23   , 'Ir' : 192.217  , 'Pt' : 195.084  ,
         'Au' : 196.9666 , 'Hg' : 200.59   , 'Tl' : 204.3833 ,
         'Pb' : 207.2    , 'Bi' : 208.9804 , 'Po' : 209.     ,
         'At' : 210.     , 'Rn' : 222.     , 'Fr' : 223.     ,
         'Ra' : 226.     , 'Ac' : 227.     , 'Th' : 232.0381 ,
         'Pa' : 231.0359 , 'U'  : 238.0289 , 'Np' : 237.     ,
         'Pu' : 244.     , 'Am' : 243.     , 'Cm' : 247.     ,
         'Bk' : 247.     , 'Cf' : 251.     , 'Es' : 252.     ,
         'Fm' : 257.     , 'Md' : 258.     , 'No' : 259.     ,
         'Lr' : 262.     , 'Rf' : 261.     , 'Db' : 262.     ,
         'Sg' : 266.     , 'Bh' : 264.     , 'Hs' : 277.     ,
         'Mt' : 268.     , 'Ds' : 281.     , 'Rg' : 272.     ,
         'Cn' : 285.     , 'Uut': 284.     , 'Uuq': 289.     ,
         'Uup': 288.     , 'Uuh': 292.     , 'Uus': 291.     ,
         'Uuo': 294.     , 'EP' : 0.000000 }

def _parse_info(form) :
    """
    Parse weight, number of atoms and number of heavy atoms from molecular formula
    """
    w = 0
    n = 0
    nh = 0
    for part in re.findall("[A-Z]+[0-9]+",form):
        m = re.match("([A-Z]+)([0-9]+)",part)
        element = m.group(1)
        number = int(m.group(2))
        w += mass[element.capitalize()]*number
        n += number
        if element != "H" : nh += number
    return w,n,nh

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to analyse the weights of solutes")
    argparser.add_argument('-db', '--db', help="the molecule database")
    argparser.add_argument('-solvent1', '--solvent1', help="the first solvent")
    argparser.add_argument('-solvent2', '--solvent2', help="the second solvent")
    argparser.add_argument('--plot',action="store_true",help="turn on plotting of distributions",default=False)
    argparser.add_argument('-cutoff','--cutoff',type=int,help="cut-off in heavy atoms",default=10)
    args = argparser.parse_args()

    db = dblib.SolvDb(filename=args.db,type="abs",filehandle="^0")

    if args.plot:
        allw = []
        alln = []
        allnh = []
    
    nfound = 0
    for entry1,entry2 in db.itersoluteoverlap(args.solvent1,args.solvent2) :
        w,n,nh = _parse_info(entry1.Formula)
        if nh > 0 and nh <= args.cutoff:
            print entry1.SoluteName+";"+entry1.Formula+";"+"%.3f"%w
            nfound += 1

        if args.plot :
            allw.append(w)
            alln.append(n)
            allnh.append(nh)

    print "____\n%d"%nfound

    if args.plot:
        for i,(data,name) in enumerate(zip([allw,alln,allnh],["allw","alln","allnh"]),1):
            f = plt.figure(i)
            h,e = np.histogram(np.asarray(data),bins=10)
            f.gca().plot((e[:-1]+e[1:])/2.0,h)
            f.savefig(name+".png",format="png")
