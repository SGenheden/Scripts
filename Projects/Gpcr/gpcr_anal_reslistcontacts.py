# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Script to calculate lifetime of specific contacts from state files
"""
import os
import sys
import argparse

import numpy as np

import pycontacts
import gpcr_lib

def _make_groups(sites, nmol, mol):

    # Load the protein template to obtain residue information
    template = gpcr_lib.load_template(mol)
    residues = np.array(template.residues)
    codes = np.array(template.codes)
    names = np.array(template.names)

    # Setup sites/Groups
    maxmem = 0
    grouplbl = []
    for site in sites:
        sidx = [template.residues.index(int(r)) for r in site.split(",")]
        maxmem = max(maxmem,len(sidx))
        grouplbl.append(",".join(["%s%d"%(codes[idx],residues[idx]) for idx in sidx]))

    groups = np.zeros((maxmem,len(sites)*nmol),dtype=int)-1
    offset = 0
    groupi = 0
    for site in sites:
        ntake = len(site.split(","))
        for mi in range(nmol):
            gidx = [offset+mi+nmol*i for i in range(ntake)]
            groups[:len(gidx),groupi] = np.asarray(gidx)
            groupi += 1
        offset += ntake*nmol

    return groups, grouplbl

def _anal_lifetimes(filename, groups, nmol, time) :
    """
    Main work routine to analyse life times

    Parameters
    ----------
    filename : string
        the state file to analyse
    groups : numpy.ndarray
        the group definitions
    nmol : int
        the number of molecules (ligands, e.g. chol) analysed
    time : float
        the total simulation time

    Returns
    -------
    list of Numpy arrays
        the results
    """

    states = gpcr_lib.read_statefile(filename)
    gstates = pycontacts.make_group_series(states,groups)

    # Conversion factor from snapshot lifetime to ns lifetime
    ns_per_snapshot = time / float(states.shape[0])
    glife_av,glife_max = pycontacts.lifetime(gstates)

    siteavs = []
    sitemaxs = []
    for igroup in range(groups.shape[1]):
        selon = glife_max[igroup*nmol:(igroup+1)*nmol] > 0
        siteavs.append([l for l in glife_av[igroup*nmol:(igroup+1)*nmol][selon]*ns_per_snapshot])
        sitemaxs.append([l for l in glife_max[igroup*nmol:(igroup+1)*nmol][selon]*ns_per_snapshot])
    return siteavs,sitemaxs

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Analysing lifetimes from state files")
    parser.add_argument('-f','--file',help="the input files.",default=[])
    parser.add_argument('--sites',nargs="+",help="sites of residues to group")
    parser.add_argument('--nmol',type=int,help="the number of molecules",default=60)
    parser.add_argument('--mol',choices=["b2","a2a","b2_a","a2a_a"],help="the protein molecules, should be either 'b2' or 'a2a'",default="b2")
    parser.add_argument('--time',type=float,help="total simulation time in ns",default=50000)
    parser.add_argument('--repeats',nargs="+",help="replacement pattern for multiple repeats",default=["r1_","r2_","r3_","r4_","r5_"])
    args = parser.parse_args()

    groups,grouplbl = _make_groups(args.sites, args.nmol, args.mol)

    siteavs,sitemaxs = _anal_lifetimes(args.file, groups, args.nmol, args.time)
    if args.repeats:
        for r in args.repeats[1:]:
            rsiteavs,rsitemaxs = _anal_lifetimes(args.file.replace(args.repeats[0],r),
                                    groups, args.nmol, args.time)
            for gi in range(groups.shape[1]):
                siteavs[gi].extend(rsiteavs[gi])
                sitemaxs[gi].extend(rsitemaxs[gi])

    for siteav,sitemax,lbl in zip(siteavs,sitemaxs,grouplbl):
        print "%s %.3f %.3f"%(lbl,np.asarray(siteav).mean(),np.asarray(sitemax).mean())
