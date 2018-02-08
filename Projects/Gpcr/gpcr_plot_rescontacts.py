# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse state files in order plot statistics of residue contacts.

It will analyse group of files, each group with a number of repeats. Three plots
will be produced
* a residue-residue contact joint probability plot, one for each group
* a residue contact probability plot, one for each group
* an amino-acid average contact probability plot, common for all groups

Examples
--------
gpcr_plot_rescontacts.py -f r1_md3_en_fit_chol.resstate.6.dat
                            r1_md3_en_fit_chol-oh.resstate.6.dat -l com oh --mol b2
"""

import os
import argparse
import sys

import numpy as np
import matplotlib
if not "DISPLAY" in os.environ or os.environ["DISPLAY"] == "" :
  matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pycontacts
import gpcr_lib
from sgenlib import colors
from sgenlib import pdb

def _plot_text(axis,text,xcoord,ycoord,arrows=False) :
  """
  Plot a series of text and try to make them non-overlapping

  Parameters
  ----------
  axis : Axis object
    the axis to do the plot on
  text : NumpyArray of string
    the text to plot
  xcoord : NumpyArray
    the original x coordinate of the texts
  ycoord : NumpyArray
    the original y coordinate of the texts
  arrows : boolean, optional
    if to draw arrows between text and original position
  """

  # Approximate text width and height
  txt_width = 0.03*(axis.get_xlim()[1] - axis.get_xlim()[0])
  txt_height = 0.085*(axis.get_ylim()[1] - axis.get_ylim()[0])

  # Find the text positions
  a = zip(ycoord, xcoord)
  text_positions = ycoord.copy()
  for index, (y, x) in enumerate(a):
    local_text_positions = [i for i in a if i[0] > (y - txt_height)
                            and (abs(i[1] - x) < txt_width * 2) and i != (y,x)]
    if local_text_positions:
      sorted_ltp = sorted(local_text_positions)
      if abs(sorted_ltp[0][0] - y) < txt_height: #True == collision
        differ = np.diff(sorted_ltp, axis=0)
        a[index] = (sorted_ltp[-1][0] + txt_height, a[index][1])
        text_positions[index] = sorted_ltp[-1][0] + txt_height
        for k, (j, m) in enumerate(differ):
          #j is the vertical distance between words
          if j > txt_height * 2: #if True then room to fit a word in
            a[index] = (sorted_ltp[k][0] + txt_height, a[index][1])
            text_positions[index] = sorted_ltp[k][0] + txt_height
            break

  # Draw the text on the modified coordinates
  bbox_props = dict(boxstyle="square", fc="w", ec="w", alpha=0.5)
  for x,y,t,tx in zip(xcoord, ycoord, text_positions,text):
    axis.text(x - txt_width, 1.01*t, tx, color='black',rotation=0,bbox=bbox_props)
    if arrows and y != t:
      axis.arrow(x, t,0,y-t, color='black',alpha=0.1, width=txt_width*0.1,
                 head_width=txt_width, head_length=txt_height*0.5,
                 zorder=0,length_includes_head=True)

def _draw_2d(axis,residues0,residues,contacts) :
  """
  Draw a residue-residue contact joint probability plot as a contour

  Parameters
  ----------
  axis : Axis object
    the axis to draw on
  residues0 : numpy array
    the tick positions of the residues
  residues : numpy array
    the x-ray number of the residues
  contacts : numpy array
    the contact probabilities
  """

  X, Y = np.meshgrid(residues0,residues0)
  axis.contour(X,Y,np.log(contacts/100.0))
  sel = residues > 0
  axis.set_xticks(residues0[sel][::20])
  axis.set_yticks(residues0[sel][::20])
  axis.set_xticklabels(residues[sel][::20])
  axis.set_yticklabels(residues[sel][::20])
  axis.set_xlabel("Residue")
  axis.set_ylabel("Residue")

def _draw_1d(axis,residues0,residues,codes,helices,contacts,ylabel,setylim,plottext=True) :
  """
  Draw a residue contact probability plot

  Parameters
  ----------
  axis : Axis object
    the axis to draw on
  residues0 : numpy array
    the tick positions of the residues
  residues : numpy array
    the x-ray number of the residues
  codes : numpy array
    the 1 amino acid residue names
  helices : list
    the first and last residue index for each helix
  contacts : numpy array
    the contact probabilities
  ylabel : string
    the label of the y-axis
  setylim : boolean, optional
    if to set ylim to  [0, 100]
  plottext : boolean, optional
    if to put out text on the 95% residues
  """

  # Plot the average residue occupancy
  axis.plot(residues0,contacts,'--k')
  for i,h in enumerate(helices) :
    c = colors.color(i) if i > 0 else colors.color(11)
    axis.plot(residues0[h[0]-1:h[1]],contacts[h[0]-1:h[1]],color=c)

  # Fix extent of axis and xticks
  axis.set_xlim([0,residues0[-1]])
  axis.set_ylim([0,setylim])
  sel = residues > 0
  axis.set_xticks(residues0[sel][::20])
  axis.set_xticklabels(residues[sel][::20])

  # Plot labels for the 95% percentile
  sel = contacts >= np.percentile(contacts,98)
  text = ["%s%d"%(c,r) for c,r in zip(codes[sel],residues[sel])]
  if plottext : _plot_text(axis,text,residues0[sel],contacts[sel])
  print "\n"+ylabel
  print "\tAverage of 95 percentile %.3f"%(contacts[sel].mean())
  print "\t95 percentile: %s"%(", ".join("%s:%.3f"%(t,c) for t,c in zip(text,contacts[sel])))

  axis.set_ylabel(ylabel)
  axis.set_xlabel("Residue")

def _draw_aa(axis,names,contacts_list,labels) :
  """
  Draw an amino-acid averaged probability plot

  Parameters
  ----------
  axis : Axis object
    the axis to draw on
  names : list of string
    the residue names
  contact_list : list of numpy array
    the residue-residue joint probability for at least one group
  labels : list of strings
    the labels of the groups
  """

  def _aa_contacts(names,contacts) :
    aas = pdb.codes.keys()
    aacount = {aa : 0 for aa in aas}
    aapp = {aa : 0 for aa in aas}
    # Accumulate probability and amino acid counts
    for name,pp in zip(names,contacts) :
      aapp[name.lower()] += pp
      aacount[name.lower()] += 1
    # Compute average probability
    for aa in aas :
      if aacount[aa] == 0 : continue
      aapp[aa] = aapp[aa] / float(aacount[aa])
    # Group histidine residues
    aapp["his"] = (aapp["his"]+aapp["hie"]+aapp["hid"]+aapp["hip"])/4.0
    del aapp["hie"]
    del aapp["hid"]
    del aapp["hip"]

    aas = sorted(aapp.keys())
    return aas,np.array([aapp[aa] for aa in aas])

  n = len(contacts_list)
  aas = []
  for i,(contacts,label) in enumerate(zip(contacts_list,labels)) :
    aas,height = _aa_contacts(names,contacts.diagonal())
    left = np.arange(1,len(aas)+1)-0.3+0.3*i
    c = colors.color(i) if i > 0 else colors.color(11)
    axis.bar(left,height,width=0.6/float(n),color=c,label=label)

  if len(labels) > 1 : axis.legend(loc='best', fancybox=True,labelspacing=0.20)
  axis.set_xticks(left+0.3)
  axis.set_xticklabels([aa.capitalize() for aa in aas],rotation='vertical')
  axis.set_xlim([0,left[-1]+1.3])
  axis.set_ylabel("Average contact probability")

def _rescontacts(filenames,labels,repeats,sites,out,mol,time,every,block=None) :
  """
  Main analysis routine

  Parameters
  ----------
  filenames : list of strings
    the group of files to analyse
  labels : list of string
    the label for each group
  repeats : list of string
    replacement pattern for multiple repeats within each group
  sites : list of string
    sites to group and analyse
  out : string
    output prefix
  mol : string
    protein identifier
  time : float
    the total simulation time
  every : int
    the reading frequency
  block : tuple of int
    just do the calculation on a block of the series
  """

  # Load the protein template to obtain residue information
  template = gpcr_lib.load_template(mol)
  residues = np.array(template.residues)
  residues0 = np.arange(1,residues.shape[0]+1)
  codes = np.array(template.codes)
  names = np.array(template.names)

  # Setup analysis of sites/groups
  if sites :
    # If the first element is l, we will do a lead-analysis
    # i.e. the first residue in a site needs to be on if the whole
    # site is to be considered to be on
    if sites[0].lower() == "l":
        makefnc = pycontacts.make_group_series_lead
        sites = sites[1:]
    else:
        makefnc = pycontacts.make_group_series

    # Find maximum member and setup group labels
    maxmem = 0
    grouplbl = []
    for site in sites:
      sidx = [template.residues.index(int(r)) for r in site.split(",")]
      maxmem = max(maxmem,len(sidx))
      grouplbl.append(",".join(["%s%d"%(codes[idx],residues[idx]) for idx in sidx]))

    # Fill the np array with group indices or -1
    groups = np.zeros((maxmem,len(sites)),dtype=int)-1
    for i,site in enumerate(sites):
        sidx = [template.residues.index(int(r)) for r in site.split(",")]
        groups[:len(sidx),i] = np.asarray(sidx)

  pcontacts = []
  pcontacts_std = []
  for fi,(filename,label) in enumerate(zip(filenames,labels)) :
    # Read the state file from disc and perform pairwise contact analysis
    states = [gpcr_lib.read_statefile(filename,every,block)]
    pp = pycontacts.pairwise_contacts(states[-1])
    if sites:
      gstates = [makefnc(states[-1],groups)]
      pg = pycontacts.pairwise_contacts(gstates[-1]).diagonal()

    # Do the same for multiple repeats and average over them
    if repeats is not None :
      pcontacts.append(np.zeros([len(repeats),pp.shape[0],pp.shape[1]]))
      pcontacts[-1][0,:,:] = pp
      if sites:
        pgroupcontacts = np.zeros([len(repeats),pg.shape[0]])
        pgroupcontacts[0,:] = pg
      for ri,r in enumerate(repeats[1:],1) :
        filename2 = filename.replace(repeats[0],r)
        states.append(gpcr_lib.read_statefile(filename2,every,block))
        pp = pycontacts.pairwise_contacts(states[-1])
        pcontacts[-1][ri,:,:] = pp
        if sites:
            gstates.append(makefnc(states[-1],groups))
            pg = pycontacts.pairwise_contacts(gstates[-1]).diagonal()
            pgroupcontacts[ri,:] = pg
      pcontacts_std.append(pcontacts[-1].std(axis=0)*100.0/np.sqrt(pcontacts[-1].shape[0]))
      pcontacts[-1] = pcontacts[-1].mean(axis=0)*100.0

      if sites:
          pgroupcontacts_std = pgroupcontacts.std(axis=0)*100.0/np.sqrt(pgroupcontacts.shape[0])
          pgroupcontacts = pgroupcontacts.mean(axis=0)*100.0
    else :
      pcontacts.append(pp*100.0)
      pgroupcontacts = pg*100

    if sites:
        print "Sites occupancies:"
        print "\n".join(["\t%s: %.3f +- %.3f"%(lbl,pg,std) for pg,std, lbl in zip(pgroupcontacts, pgroupcontacts_std, grouplbl)])

    # Conversion factor from snapshot lifetime to ns lifetime
    print time
    if block is not None:
        time = time / float(block[0]) * block[1]
        print time
    ns_per_snapshot = time / float(states[0].shape[0])

    # Calculate lifetimes
    all_states = np.concatenate(states,axis=1)
    life_av,life_max = pycontacts.lifetime(all_states)
    life_av = np.reshape(life_av,[len(repeats),residues.shape[0]]).mean(axis=0)*ns_per_snapshot
    life_max = np.reshape(life_max,[len(repeats),residues.shape[0]]).mean(axis=0)*ns_per_snapshot
    print "Median av-lifetime = %.3f"%np.median(life_av)
    print "Average av-lifetime = %.3f"%np.mean(life_av)
    print "Median max-lifetime = %.3f"%np.median(life_max)
    print "Average max-lifetime = %.3f"%np.mean(life_max)

    if sites:
        all_gstates = np.concatenate(gstates,axis=1)
        glife_av,glife_max = pycontacts.lifetime(all_gstates)
        glife_av = np.reshape(glife_av,[len(repeats),groups.shape[1]]).mean(axis=0)*ns_per_snapshot
        glife_max = np.reshape(glife_max,[len(repeats),groups.shape[1]]).mean(axis=0)*ns_per_snapshot
        print "Sites occ. av-lifetime max-lifetime:"
        print "\n".join(["\t%s: %.3f\t%.3f\t%.3f"%(lbl,o,a,m) for o,a,m,lbl in zip(pgroupcontacts,glife_av,glife_max,grouplbl)])

    # Draw a 2D residue-residue joint probability plot
    f2d = plt.figure(10+fi)
    _draw_2d(f2d.gca(),residues0,residues,pcontacts[-1])
    f2d.savefig("%s_%s_2d.png"%(out,label),format="png",dpi=300)

    # Draw a residue contact probability plot
    f1d = plt.figure(20+fi,figsize=(6.85,3.41),tight_layout=True)
    p = 80 if label == "com" else 100
    _draw_1d(f1d.gca(),residues0,residues,codes,template.rhelices,pcontacts[-1].diagonal(),"Contact probability",p)
    f1d.savefig("%s_%s_1d.png"%(out,label),format="png",dpi=300)

    # Draw a average lifetime plot
    f1d = plt.figure(40+fi,figsize=(6.85,3.41),tight_layout=True)
    _draw_1d(f1d.gca(),residues0,residues,codes,template.rhelices,life_av,"Average lifetime (ns)",0.6)
    f1d.savefig("%s_%s_avlife.png"%(out,label),format="png",dpi=300)

    # Draw a maximum lifetime plot
    f1d = plt.figure(50+fi,figsize=(6.85,3.41),tight_layout=True)
    _draw_1d(f1d.gca(),residues0,residues,codes,template.rhelices,life_max,"Maximum lifetime (ns)",35)
    f1d.savefig("%s_%s_maxlife.png"%(out,label),format="png",dpi=300)

    # And print it out do disc
    with open("%s_%s_1d.txt"%(out,label),"w") as f :
      for (name,res,prob,prob_std,rav,rmax) in zip(names,residues,pcontacts[-1].diagonal(),pcontacts_std[-1].diagonal(),life_av,life_max) :
        f.write("%s%d\t%8.3f\t%8.3f\t%8.3f\t%8.3f\n"%(name.capitalize(),res,prob,prob_std,rav,rmax))

  # Plot residue-type averaged occupancies
  faa = plt.figure(30)
  _draw_aa(faa.gca(),names,pcontacts,labels)
  faa.savefig(args.out+"_aa.png",format="png",dpi=300)

if __name__ == '__main__' :

  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Plotting residue contacts")
  parser.add_argument('-f','--files',nargs='+',help="a list of input files.",default=[])
  parser.add_argument('-l','--labels',nargs='+',help="a label for each input file.",default=[])
  parser.add_argument('-o','--out',help="the output prefix.",default="rescontacts")
  parser.add_argument('--mol',choices=["b2","a2a","b2_a","a2a_a"],help="the protein molecules, should be either 'b2' or 'a2a'",default="b2")
  parser.add_argument('--repeats',nargs="+",help="replacement pattern for multiple repeats",default=["r1_","r2_","r3_","r4_","r5_"])
  parser.add_argument('--sites',nargs="+",help="sites of residues to group")
  parser.add_argument('--time',type=float,help="total simulation time in ns",default=50000)
  parser.add_argument('--every',type=int,help="reading frequency",default=1)
  parser.add_argument('--block',type=int,nargs="+",help="block analysis")
  args = parser.parse_args()

  _rescontacts(args.files,args.labels,args.repeats,args.sites,args.out,args.mol,args.time,args.every,args.block)
