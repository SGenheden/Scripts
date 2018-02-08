# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse CG simulations of GPCRs.

It can perform the following kind of analysis
* buried   - analyse how many cholesterol is found when the cut-off for buriedeness
             is varied.
* contacts - analyse cholesterol in contacts with the protein. The cut-off
             is set with the --cutoff argument. Writes out three states files
             a. cholesterol in contact with the protein, b. cholesterol in contact
             with all residues in the protein, c. buried cholesterol in contact
             with all residues in the protein. These states files are written
             out both as measured from the centroid of the cholesterol AND
             as measured from the OH group.
* lipidcontacts - similar to contacts but for the short and long tails of POPC.
                  independent analysis carried out for each chain. Does not
                  look for buried lipids.
* counting - count the number of cholesterol/lipids in the two leaflets and
             number of buried cholesterols. Write out these to state files.
* density  - discretises the cholesterol/lipids in the membrane plane on a 2D grid.
             Create densities for a. lipid head groups (head) b. lipid centroid (popc)
             c. glycerol beads (glycerol) d. select tail beads (tails) e. cholesterol
             centroid (chol) f. cholesterol excluding buried (chol2) g. hydroxyl bead
             of cholesterol (hydroxyl) h. same as g, but excluding buried (hydroxyl2)
             i. centroid of buried cholesterol (mid_chol2)
             j. hydroxyl of buried cholesterol (mid_hydroxyl2)
             for the densities a.-g. a density for the lower and upper leaflet is created,
             prefixed with low_ and upp_. In addition to the density this analysis
             also keeps track of the total amount of discretised molecules in each
             grid, the number of snapshots and the box size. All information is saved
             to an .npz-file. The density at specific time point during the simulation
             is also written out.
* order  - analyse the order parameters of the lipids. Writes out information to an
           .npz-file
* thickness - analyse the bilayer thickness and writes out information to an .npz-file
* savemol - print out cholesterols and lipids in contacts with the protein to xyz-files

Examples
--------
gpcr_mdanal.py -x r1_md3.xtc -s r1_md2.gro -a counting density contacts
gpcr_mdanal.py -x r1_md3.xtc -s r1_md2.gro -p out/ -a density
gpcr_mdanal.py -x r1_md3.xtc -s r1_md2.gro -c 5.0
"""

import argparse
import ConfigParser
import sys
import os

import numpy as np
import MDAnalysis as md
import MDAnalysis.analysis as anal
import MDAnalysis.lib.util as mdcore
md.core.flags['use_periodic_selections'] = True
md.core.flags['use_KDTree_routines'] = False

import gpcr_lib

def _get_selections(filename) :

  selconfig = {}
  if filename is not None :
    p = ConfigParser.SafeConfigParser()
    p.read(filename)
    for o in p.options("Selections") :
      selconfig[o] = p.get("Selections",o)
  else :
    selconfig["protein"]="protein"
    selconfig["phosphate"]="name PO4"
    selconfig["choline"]="name NC3"
    selconfig["hydroxyl"]="name ROH"
    selconfig["glycerol1"]="name GL2"
    selconfig["glycerol2"]="name GL1"
    #   This were the name in old MARTINI POPC
    #selconfig["short"]="name C2A,name C3A"
    #selconfig["long"]="name C2B,name D3B,name C4B"
    selconfig["short"]="name C2B,name C3B"
    selconfig["long"]="name D2A,name C3A,name C4A"
    selconfig["lipid"]="resname POPC"
    selconfig["cholesterol"]="resname CHOL"
    selconfig["lipid-around"]="byres ((name C2B or name C3B or name D2A or name C3A or name C4A) and around %.3f protein)"
    selconfig["chol-around"]="resname CHOL and around %.3f protein"
    selconfig["hydroxyl-around"]="byres (name ROH and around %.3f protein)"
    selconfig["bond-names"]="C1B-C2B C2B-C3B C3B-C4B C1A-D2A D2A-C3A C3A-C4A"
  return selconfig

if __name__ == '__main__' :

  # Command-line input
  parser = argparse.ArgumentParser(description="Analyzing GPCR simualtions")
  parser.add_argument('-x','--xtc',help="the xtc file.",default="")
  parser.add_argument('-s','--struct',help="a structure file",default="r1_md2_fit.gro")
  parser.add_argument('-p','--prefix',help="an output prefix",default="")
  parser.add_argument('-a','--analysis',nargs="+",choices=["density","counting","joint","thickness","order","savemol","contacts","lipidcontacts","buried"],help="what analysis do to",default=[])
  parser.add_argument('-t','--time',type=int,help="ps per snapshot",default=20)
  parser.add_argument('-b','--box',type=int,choices=[52,35],help="the truncated box size",default=35)
  parser.add_argument('-c','--cutoff',type=float,help="contact cut-off",default=6.0)
  parser.add_argument('-bc','--bcutoff',type=float,help="buried cut-off",default=8.0)
  parser.add_argument('--skip',type=int,help="skip every")
  parser.add_argument('--selections',help="a configuration file with custom selections")
  parser.add_argument('--reslist',help="a list of residues to analyse individual chol contacts")
  args = parser.parse_args()

  do_density = "density" in args.analysis
  do_counting = "counting" in args.analysis
  do_thickness = "thickness" in args.analysis
  do_order = "order" in args.analysis
  do_savemol = "savemol" in args.analysis
  do_contacts = "contacts" in args.analysis
  do_lipidcontacts = "lipidcontacts" in args.analysis
  do_countburied = "buried" in args.analysis
  do_joint = "joint" in args.analysis
  print "Doing %s analysis"%", ".join(args.analysis)

  # Make output filenames
  h,t = os.path.splitext(args.xtc)
  prefix = args.prefix+h

  # Create an MDAnalysis universe
  universe=md.Universe(args.struct,args.xtc)

  # Make selections of interest
  selconfig = _get_selections(args.selections)
  prot = universe.select_atoms(selconfig["protein"])
  residues = prot.residues
  heads = universe.select_atoms(selconfig["phosphate"]) # Lipids head group
  hydroxyl = universe.select_atoms(selconfig["hydroxyl"]) # Cholesterol hydroxyl bead
  glycerols = universe.select_atoms(selconfig["glycerol1"]) # Glycerol group
  glycerols2 = universe.select_atoms(selconfig["glycerol2"]) # Glycerol group1
  nas = universe.select_atoms(selconfig["choline"]) # Lipids head group
  short_chain = [universe.select_atoms(n) for n in selconfig["short"].split(",")]
  long_chain = [universe.select_atoms(n) for n in selconfig["long"].split(",")]
  lipids= [res for res in universe.select_atoms(selconfig["lipid"]).residues] # List of all lipids
  chols = [res for res in universe.select_atoms(selconfig["cholesterol"]).residues] # List of all cholesterols
  if do_savemol :
    lipid_around = universe.select_atoms(selconfig["lipid-around"]%args.cutoff)
    if len(chols) > 0 :
      chol_around = universe.select_atoms(selconfig["chol-around"]%args.cutoff)
      hydroxyl_around = universe.select_atoms(selconfig["hydroxyl-around"]%args.cutoff)
  # Lipids bond vectors (the ones below are for POPC)
  bond_sel = []
  for (i,b) in enumerate(selconfig["bond-names"].split()) :
    start,end = b.split("-")
    bond_sel.append([universe.select_atoms("name %s"%start),universe.select_atoms("name %s"%end)])

  # Obtain the first atom index of each residue, first index is 0
  rfirst = np.zeros(len(residues)+1,dtype=np.int)
  rfirst[i] = 0
  for i,residue in enumerate(residues[1:],1) :
    rfirst[i] = residue[0].number-residues[0][0].number
  rfirst[-1] = prot[-1].number-residues[0][0].number

  # Make grid data
  grid = gpcr_lib.GridTemplate(universe.coord._pos)

  # Allocated count matrices
  if do_counting :
    lipleaffile = open(prefix+"_lip.leaflet.dat","wb")
    if len(chols) > 0 :
      cholcount_low = 0.0
      cholcount_upp = 0.0
      cholcount_low2 = 0.0
      cholcount_upp2 = 0.0
      cholmidfile = open(prefix+"_chol.mid.dat","wb")
      cholleaffile = open(prefix+"_chol.leaflet.dat","wb")

  # Allocate data for density analysis
  if do_density :
    density_mat = {}
    density_stride = 500000  / args.time  #10000000
    density_names = "head popc glycerol tails chol chol2 hydroxyl hydroxyl2".split()
    for name in density_names :
      density_mat["low_"+name] = grid.create()
      density_mat["upp_"+name] = grid.create()
    density_mat["mid_hydroxyl2"] = grid.create()
    density_mat["mid_chol2"] = grid.create()
    density_mat["cnt_popc"] = [0,0]
    density_mat["cnt_chol"] = [0,0]
    density_mat["cnt_chol2"] = [0,0,0]
    density_mat["nsnap"] = 0
    density_mat["box"] = 2*args.box

  # Allocate data for thickness analysis
  if do_thickness :
    thickness_mat = {}
    thickness_mat["lower"] = grid.create()+1000
    thickness_mat["upper"] = grid.create()+1000

  # Allocate data for order parameter analysis
  if do_order :
    order_mat = {}
    order_mat["all_sum"] = np.zeros([len(bond_sel),3])
    order_mat["low_sum"] = grid.create()
    order_mat["low_count"] = grid.create()
    order_mat["upp_sum"] = grid.create()
    order_mat["upp_count"] = grid.create()

  # Open files for savemol analysis
  if do_savemol :
    save_stride = 25000 / args.time
    crd_lip = {"file":open(prefix+"_lip.xyz","w"),"frame":0}
    if len(chols) > 0 :
      crd_chol = {"file":open(prefix+"_chol.xyz","w"),"frame":0}
      crd_chol2 = {"file":open(prefix+"_chol_hydroxyl.xyz","w"),"frame":0}

  # Open files for contact analysis
  if do_lipidcontacts :
    shortmolfile = open(prefix+"_short.mstate.%.0f.dat"%args.cutoff,"wb")
    shortresfile = open(prefix+"_short.resstate.%.0f.dat"%args.cutoff,"wb")
    longmolfile = open(prefix+"_long.mstate.%.0f.dat"%args.cutoff,"wb")
    longresfile = open(prefix+"_long.resstate.%.0f.dat"%args.cutoff,"wb")

  if do_contacts and len(chols) > 0 :
    cholmolfile = open(prefix+"_chol.mstate.%.0f.dat"%(args.cutoff),"wb")
    cholresfile = open(prefix+"_chol.resstate.%.0f.dat"%(args.cutoff),"wb")
    cholburresfile = open(prefix+"_chol.buried.rstate.%.0f.dat"%(args.cutoff),"wb")
    ohmolfile = open(prefix+"_chol-oh.mstate.%.0f.dat"%(args.cutoff),"wb")
    ohresfile = open(prefix+"_chol-oh.resstate.%.0f.dat"%(args.cutoff),"wb")
    ohburresfile = open(prefix+"_chol-oh.buried.rstate.%.0f.dat"%(args.cutoff),"wb")
    if args.reslist is not None:
      with open(args.reslist,"r") as f :
        lines = f.readlines()
      template = gpcr_lib.load_template(lines[0].strip())
      reslist = [template.residues.index(int(l.strip())) for l in lines[1:]]
      reslistfile = open(prefix+"_chol-oh.resliststate.%.0f.dat"%(args.cutoff),"wb")
    else :
      reslist = None
      reslistfile = None


  if do_joint and len(chols) > 0 :
    jointcom = np.zeros([len(residues),len(residues)])
    jointcomburr = np.zeros([len(residues),len(residues)])
    jointoh = np.zeros([len(residues),len(residues)])
    jointohburr = np.zeros([len(residues),len(residues)])

  # Allocate matrices for buried analysis
  if do_countburied :
    buried_depth = [2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,12.0]
    nburied = np.zeros(len(buried_depth))

  #
  # Main loop
  #
  nfiles=1
  nsnap = 0
  for ti,ts in enumerate(universe.trajectory) :
    if args.skip is not None and (ti % args.skip == 0) : continue
    nsnap += 1

    # Calculates some properties important for many analysis
    box = ts.dimensions[:3]
    # This is taken as a the middle of the bilayer
    midz = heads.get_positions()[:,2].mean()
    # Discretise head group positions
    head_idx   = grid.indices(heads.get_positions())
    # Find center of lipids, glycerol beads, and tail chains
    lipid_cent = np.array([lipid.centroid() for lipid in lipids])
    glyc_cent = (glycerols.get_positions() + glycerols2.get_positions())/2.0
    short_chain_cent = (short_chain[0].get_positions() + short_chain[1].get_positions())/2.0
    long_chain_cent = (long_chain[0].get_positions() + long_chain[1].get_positions() + long_chain[2].get_positions())/3.0
    both_chain_cent = (short_chain[0].get_positions() + short_chain[1].get_positions()+long_chain[0].get_positions() + long_chain[1].get_positions() + long_chain[2].get_positions())/5.0

    if do_counting : gpcr_lib.write_booleans(lipleaffile,heads.get_positions()[:,2]<midz)
    if len(chols) > 0 :
      # Calculates cholesterol centers
      chol_cent  = np.array([chol.centroid() for chol in chols])
      # Find and write out buried cholesterols
      buried = np.abs(hydroxyl.get_positions()[:,2]-midz)<=args.bcutoff
      notburied = np.logical_not(buried)
      if do_counting :
        gpcr_lib.write_booleans(cholmidfile,buried)
        gpcr_lib.write_booleans(cholleaffile,chol_cent[:,2]<midz)
        # Count number of cholesterols
        cholcount_low = cholcount_low+(chol_cent[:,2]<midz).sum()
        cholcount_upp = cholcount_upp+(chol_cent[:,2]>=midz).sum()
        cholcount_low2 = cholcount_low2+np.logical_and(chol_cent[:,2]<midz,notburied).sum()
        cholcount_upp2 = cholcount_upp2+np.logical_and(chol_cent[:,2]>=midz,notburied).sum()

      if do_countburied :
        for i,depth in enumerate(buried_depth)  :
          nburied[i] = nburied[i] + (np.abs(hydroxyl.get_positions()[:,2]-midz)<=depth).sum()
      if do_contacts :
        gpcr_lib.anal_contacts(prot,rfirst,chol_cent,args.cutoff,midz,box,cholmolfile,cholresfile,cholburresfile,buried)
        gpcr_lib.anal_contacts(prot,rfirst,hydroxyl.get_positions(),args.cutoff,midz,box,ohmolfile,ohresfile,ohburresfile,buried,reslist,reslistfile)
      if do_joint :
        jointcom,jointcomburr = gpcr_lib.anal_jointdist(prot,rfirst,chol_cent,args.cutoff,midz,box,jointcom,jointcomburr,buried)
        jointoh,jointohburr = gpcr_lib.anal_jointdist(prot,rfirst,hydroxyl.get_positions(),args.cutoff,midz,box,jointoh,jointohburr,buried)
    else :
      chol_cent = None
      buried = None

    if do_density :
      density_mat["nsnap"] = density_mat["nsnap"] + 1
      gpcr_lib.anal_density(low_sel=heads.get_positions()[:,2] < midz,count="popc",grid=grid,mat=density_mat,head=head_idx,popc=lipid_cent,glycerol=glyc_cent,tails=both_chain_cent)
      if len(chols) > 0 :
        gpcr_lib.anal_density(low_sel=chol_cent[:,2]<midz,count="chol",grid=grid,mat=density_mat,chol=chol_cent,hydroxyl=hydroxyl.get_positions())
        gpcr_lib.anal_density(low_sel=chol_cent[:,2]<midz,count="chol2",grid=grid,mat=density_mat,buried=buried,chol2=chol_cent,hydroxyl2=hydroxyl.get_positions())

      if ti %  density_stride == 0  : #
        density_file = "%s_densities%d.npz"%(prefix,nfiles)
        grid.save(density_mat,args.box,density_file)
        nfiles = nfiles + 1

    if do_order :
      gpcr_lib.anal_orderparams(heads,head_idx,midz,bond_sel,order_mat)

    if do_thickness :
      gpcr_lib.anal_thickness(heads,head_idx,midz,thickness_mat)

    if do_savemol and ti % save_stride == 0 :
      gpcr_lib.anal_savemols(lipid_around,ti,crd_lip)
      if len(chols) > 0 :
        gpcr_lib.anal_savemols(chol_around,ti,crd_chol)
        gpcr_lib.anal_savemols(hydroxyl_around,ti,crd_chol2)

    if do_lipidcontacts :
      contact_short = gpcr_lib.anal_contacts(prot,rfirst,short_chain_cent,args.cutoff,midz,box,shortmolfile,shortresfile)
      contact_long = gpcr_lib.anal_contacts(prot,rfirst,long_chain_cent,args.cutoff,midz,box,longmolfile,longresfile)

    if ti % 10001 == 0:
      print "%d"%ti,
      sys.stdout.flush()

  # Finalisation
  # Deallocation, closing of files, saving of grids

  if do_density :
    grid.save(density_mat,args.box,prefix+"_densities.npz")

  if do_thickness :
    thickness_mat["lower"][mat["lower"]==1000] = 0
    thickness_mat["upper"][mat["upper"]==1000] = 0
    grid.save(thickness_mat,args.box,prefix+"_thickness.npz")

  if do_order :
    grid.save(order_mat,args.box,prefix+"_order.npz")

  if do_savemol :
    crd_lip["file"].close()
    if len(chols) > 0 :
      crd_chol["file"].close()
      crd_chol2["file"].close()

  if do_contacts and len(chols) > 0:
    cholmolfile.close()
    cholresfile.close()
    cholburresfile.close()
    ohmolfile.close()
    ohresfile.close()
    ohburresfile.close()

  if do_lipidcontacts :
    shortmolfile.close()
    shortresfile.close()
    longmolfile.close()
    longresfile.close()

  if do_joint :
    jointcom = jointcom / float(nsnap)
    jointcomburr = jointcomburr / float(nsnap)
    jointoh = jointoh / float(nsnap)
    jointohburr = jointohburr / float(nsnap)
    np.savez(prefix+"_joint.%.0f.npz"%args.cutoff,jointcom=jointcom,jointcomburr=jointcomburr,jointoh=jointoh,jointohburr=jointohburr)

  if do_counting :
    lipleaffile.close()
    if len(chols) > 0 :
      cholmidfile.close()
      cholleaffile.close()
      print "\nAverage number of cholesterols, lower=%.3f, upper=%.3f"%(cholcount_low/float(nsnap+1),cholcount_upp/float(nsnap+1))
      print "Average number of cholesterols (excluding buried), lower=%.3f, upper=%.3f"%(cholcount_low2/float(nsnap+1),cholcount_upp2/float(nsnap+1))
      if do_countburied :
          for i,depth in enumerate(buried_depth)  :
            print "Average number of cholesterols at depth %.2f = %.3f"%(depth,nburied[i]/float(nsnap+1))

  print nsnap
