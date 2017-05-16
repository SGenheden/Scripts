
#
# THIS CODE DOES NOT WORK DUE TO PRECISION
# USE THE C++ code instead
#

import argparse

import numpy as np

from sgenlib import parsing

def _make_coulomb_ff(r, dielectric) :


    fx = np.power(r, -1.0)
    fdx = np.power(np.power(r,2.0), -1.0)

    return fx / dielectric, fdx / dielectric

def _make_rf_ff(r, cutoff, dielectric, epsilon) :

    cutoff2 = cutoff*cutoff
    cutoff3 = cutoff2*cutoff
    krf = (epsilon - dielectric) / (cutoff3*(2.0*epsilon + dielectric))
    crf = (3.0*epsilon) / (cutoff * (2.0*epsilon + dielectric))
    fshift = (3.0*dielectric) / (cutoff2 * (2.0*epsilon + dielectric))

    rinv = 1.0 / r
    rinv2 = rinv*rinv
    r2 = r * r

    fx = rinv + krf*r2 - crf
    fdx = rinv2 - 2.0*krf*r - fshift

    return fx / dielectric, fdx / dielectric

def _make_lj_ff(r) :

    r = np.power(r,-1.0, dtype=np.float64)
    r2 = np.power(r, 2.0, dtype=np.float64)
    r6 = np.power(r2, 3.0, dtype=np.float64)
    r12 = np.power(r6, 2.0, dtype=np.float64)
    r7 = np.multiply(r6, r, dtype=np.float64)
    r13 = np.multiply(r12, r, dtype=np.float64)

    gx = -r6
    gdx = np.multiply(r7,-6.0)

    hx = r12
    hdx = np.multiply(r13, 12.0)

    return gx, gdx, hx, hdx

def _make_shifted_ff(r, cutoff) :

    shiftlj = 0.9
    rmodlj = cutoff - shiftlj
    rmodlj2 = rmodlj*rmodlj
    rmodlj3 = rmodlj2*rmodlj
    rmodlj4 = rmodlj2*rmodlj2

    a2 = -(10.0*cutoff - 7.0*shiftlj)/ ((cutoff**8.0)*rmodlj2)
    b2 = (9.0*cutoff - 7.0*shiftlj)/ ((cutoff**8.0)*rmodlj3)
    c2 = 1.0/(cutoff**6.0) - 2.0*a2*rmodlj3 - 2.0*b2*rmodlj4
    a3 = -(16.0*cutoff - 13.0*shiftlj)/ ((cutoff**14.0)*rmodlj2)
    b3 = (15.0*cutoff - 13.0*shiftlj)/ ((cutoff**14.0)*rmodlj3)
    c3 = 1.0/(cutoff**12.0) - 4.0*a3*rmodlj3 - 4.0*b3*rmodlj4

    rlj = r - shiftlj
    rlj2 = rlj*rlj
    rlj3 = rlj2*rlj
    rlj4 = rlj2*rlj2


    gx = -1.0/(r**6.0) + c2
    gdx = -6.0/(r**7.0)
    hx = 1.0/(r**12.0) - c3
    hdx = 12.0/(r**13.0)

    sel = r > shiftlj
    gx[sel] = -1.0/(r[sel]**6.0) + 2.0*a2*rlj3[sel] + 1.5*b2*rlj4[sel] + c2
    gdx[sel] = -6.0/(r[sel]**7.0) - 6.0*(a2*rlj2[sel] + b2*rlj3[sel])
    hx[sel] = 1.0/(r[sel]**12.0)  - 4.0*a3*rlj3[sel] - 3.0*b3*rlj4[sel] - c3
    hdx[sel] = 7.0/(r[sel]**13.0) + 12.0*(a3*rlj2[sel] + b3*rlj3[sel])

    return gx, gdx, hx, hdx

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Script to make tabulated potentials")
    argparser.add_argument('-e','--electrostatics',choices=["coulomb","rf"],help="the electrostatic force field",default="coulomb")
    argparser.add_argument('-v','--vanderwaals',choices=["lj","shifted"],help="the van der Waals force field",default="lj")
    argparser.add_argument('--cutoff',type=float,help="the cut-off",default=1.0)
    argparser.add_argument('--dielectric',type=float,help="the dielectric constant",default=1.0)
    argparser.add_argument('--epsilon',type=float,help="the dielectric constant",default=1.0)
    argparser.add_argument('-o', '--out', help="the output filename", default="table.xvg")
    args = argparser.parse_args()

    if args.electrostatics == "coulomb" and args.vanderwaals == "lj" :
        table = parsing.parse2ndarray("/Users/samuel/Documents/Programs/Gromacs/5.1/share/gromacs/top/table6-12.xvg")
    elif args.electrostatics == "rf" and args.vanderwaals == "shifted" :
        table = parsing.parse2ndarray("/Users/samuel/Downloads/hybrid-tutorial/OPLS-MARTINI/FORCEFIELDS/mar-table_CG_CG.xvg")

    delta = 0.002000
    r = np.arange(delta, args.cutoff+2.0*delta+1.0, delta, dtype=np.float64)
    #r = np.array(table[:,0], dtype=np.float64)
    #table[r>args.cutoff,1:] = 0.0
    #r = r[1:]

    if args.electrostatics == "coulomb" :
        fx, fdx = _make_coulomb_ff(r, args.dielectric)
        fx[r<4.0E-2] = 0.0
        fdx[r<4.0E-2] = 0.0
    elif args.electrostatics == "rf" :
        fx, fdx = _make_rf_ff(r, args.cutoff, args.dielectric, args.epsilon)

    if args.vanderwaals == "lj" :
        gx, gdx, hx, hdx = _make_lj_ff(r)
        gx[r<4.0E-2] = 0.0
        gdx[r<4.0E-2] = 0.0
        hx[r<4.0E-2] = 0.0
        hdx[r<4.0E-2] = 0.0
    if args.vanderwaals == "shifted" :
        gx, gdx, hx, hdx = _make_shifted_ff(r, args.cutoff)
    """
    print "fx: ", np.max(table[1:,1]-fx),np.argmax(table[1:,1]-fx),r[np.argmax(table[1:,1]-fx)]
    print "fdx: ", np.max(table[1:,2]-fdx),np.argmax(table[1:,2]-fdx),r[np.argmax(table[1:,2]-fdx)]
    print "gx: ", np.max(table[1:,3]-gx),np.argmax(table[1:,3]-gx),r[np.argmax(table[1:,3]-gx)]
    print "gdx: ", np.max(table[1:,4]-gdx),np.argmax(table[1:,4]-gdx),r[np.argmax(table[1:,4]-gdx)]
    print "hx: ", np.max(table[1:,5]-hx),np.argmax(table[1:,5]-hx),r[np.argmax(table[1:,5]-hx)]
    print "hdx: ", np.max(table[1:,6]-hdx),np.argmax(table[1:,6]-hdx),r[np.argmax(table[1:,6]-hdx)]
    """
    """print "Coul"
    print np.abs(table[121:131,1]-fx[120:130])/np.abs(table[121:131,1])*100
    print ""
    print np.abs(table[21:31,1]-fx[20:30])/np.abs(table[21:31,1])*100
    print ""
    print np.abs(table[121:131,2]-fdx[120:130])/np.abs(table[121:131,2])*100
    print ""
    print np.abs(table[21:31,2]-fdx[20:30])/np.abs(table[21:31,2])*100

    print "r6"
    print np.abs(table[121:131,3]-gx[120:130])/np.abs(table[121:131,3])*100
    print ""
    print np.abs(table[21:31,3]-gx[20:30])/np.abs(table[21:31,3])*100
    print ""
    print np.abs(table[121:131,4]-gdx[120:130])/np.abs(table[121:131,4])*100
    print ""
    print np.abs(table[21:31,4]-gdx[20:30])/np.abs(table[21:31,4])*100

    print "r12"
    print np.abs(table[121:131,5]-hx[120:130])/np.abs(table[121:131,5])*100
    print ""
    print np.abs(table[21:31,5]-hx[20:30])/np.abs(table[21:31,5])*100
    print ""
    print np.abs(table[121:131,6]-hdx[120:130])/np.abs(table[121:131,6])*100
    print ""
    print np.abs(table[21:31,6]-hdx[20:30])/np.abs(table[21:31,6])*100"""

    with open(args.out, "w") as f :
        f.write("%20.10E %20.10E %20.10E %20.10E %20.10E %20.10E %20.10E\n"%(0, 0, 0, 0 ,0, 0, 0))
        for ri, fxi, fdxi, gxi, gdxi, hxi, hdxi  in zip(r, fx, fdx, gx, gdx, hx, hdx) :
            if ri <= args.cutoff :
                f.write("%20.10E %20.10E %20.10E %20.10E %20.10E %20.10E %20.10E\n"
                        %(ri, fxi, fdxi, gxi, gdxi, hxi, hdxi))
            else :
                f.write("%20.10E %20.10E %20.10E %20.10E %20.10E %20.10E %20.10E\n"
                        %(ri, 0, 0, 0 ,0, 0, 0))
