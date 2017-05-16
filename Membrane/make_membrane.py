import sys
import pdb
import argparse
import geo,fitting
import numpy as np
import numpy.random as random

def rotate_lipid(lipid) :

  lipidmasses = lipid.collect("masses")
  xyz = lipid.collect("xyz")
  center = xyz.mean(axis=0)
  xyz = xyz - center
  moi = geo.momentOfInertia(xyz,lipidmasses)
  princip = geo.principalAxes(moi)
  normv = np.array([0.0,0.0,1.0])
  rotvec = geo.rotaxis(princip[0,:],normv)
  alpha = geo.angle(princip[0,:],normv)
  rotmat = geo.rotation_matrix(alpha,rotvec)
  lipid.update_xyz(fitting.rotate(xyz,rotmat)+center)

def origin_lipid(lipid) :
  
  xyz = lipid.collect("xyz")
  center = xyz.mean(axis=0)
  xyz = xyz - center
  xyz[:,2] = xyz[:,2] - (xyz[:,2].max()-xyz[:,2].min())/2.0
  lipid.update_xyz(xyz)

def create_partner(lipid) :

  new = pdb.Residue()
  new.copy(lipid)
  xyz = new.collect("xyz")
  xyz[:,2] = -xyz[:,2]
  new.update_xyz(xyz)
  return new
  

if __name__ == "__main__":


  # Setup a parser of the command-line arguments
  parser = argparse.ArgumentParser(description="Program to make an initial membrane configuration")
  parser.add_argument('-l','--lipid',help="the lipid input file")
  parser.add_argument('-w','--water',help="the water input file")
  parser.add_argument('-o','--out',help="the output filename, default='membrane.pdb'",default="membrane.pdb")
  args = parser.parse_args()

  lipid = pdb.PDBFile()
  lipid.read(args.lipid,gro=args.lipid[-3:].lower()=="gro")
  
  rotate_lipid(lipid.residues[0]) 
  origin_lipid(lipid.residues[0])
  lipid.residues.append(create_partner(lipid.residues[0]))

  xyz = lipid.residues[0].collect("xyz")
  xwidth = xyz[:,0].max()-xyz[:,0].min()
  ywidth = xyz[:,1].max()-xyz[:,1].min()
  width = max(xwidth,ywidth)

  membrane = pdb.PDBFile()
  nres = 0
  for xi in range(0,8) :
    for yi in range(0,8) : 
      nres = nres + 2
      xdel = xi*xwidth      
      ydel = yi*ywidth
      new1 = pdb.Residue()
      new1.copy(lipid.residues[0])
      xyz = new1.collect("xyz")
      xyz[:,0] = xyz[:,0] + xdel + random.randn()
      xyz[:,1] = xyz[:,1] + ydel + random.randn()
      xyz[:,2] = xyz[:,2] + random.randn()
      cent = xyz.mean(axis=0)
      xyz = xyz - cent
      rotmat = geo.rotation_matrix(random.randin(0,15)*np.pi/180.0,np.array([0.0,0.0,1.0]))
      xyz = fitting.rotate(xyz,rotmat)+cent
      new1.update_xyz(xyz)
      new1.serial = nres - 1
      new2 = pdb.Residue()
      new2.copy(lipid.residues[1])
      new2.serial = nres
      xyz = new2.collect("xyz")
      xyz[:,0] = xyz[:,0] + xdel + random.randn()
      xyz[:,1] = xyz[:,1] + ydel + random.randn()
      xyz[:,2] = xyz[:,2] + random.randn()
      cent = xyz.mean(axis=0)
      xyz = xyz - cent
      rotmat = geo.rotation_matrix(random.randint(0,15)*np.pi/180.0,np.array([0.0,0.0,1.0]))
      xyz = fitting.rotate(xyz,rotmat)+cent
      new2.update_xyz(xyz)
      membrane.residues.append(new1)
      membrane.residues.append(new2)
  natom = 0
  memxyz = np.zeros([nres*len(lipid.residues[0].atoms),3])
  for residue in membrane.residues :
    for atom in residue.atoms :
      natom = natom + 1
      atom.serial = natom
      atom.residue = residue.serial
      memxyz[natom-1,:] = atom.xyz
  membrane.write(args.out)
  mem_lenx = memxyz[:,0].max() - memxyz[:,0].min()
  mem_leny = memxyz[:,1].max() - memxyz[:,1].min()
  
  watbox = pdb.PDBFile(filename=args.water)
  
  nwattot = 40*nres
  nwatbox = len(watbox.residues)
  nbox    = int(np.ceil(float(nwattot)/float(nwatbox)))
  toadd   = nbox*nwatbox/2
  
  watdel = delta  = (1.0/0.0335)**0.3333333
  print watdel
  watdel2 = watdel / 2.0
  x = memxyz[:,0].min()
  y = memxyz[:,1].min()
  z1 = memxyz[:,2].max()-watdel2
  z2 = memxyz[:,2].min()+watdel2
  for i in range(toadd) :
    nres = nres + 1
    idx,idx2 = random.randint(0,len(watbox.residues),size=2)
    new1 = pdb.Residue()
    new1.copy(watbox.residues[idx])
    xyz = new1.collect("xyz")
    xyz = xyz - xyz[0,:]
    xyz[:,0] = xyz[:,0] + x + watdel2 + 0.2*random.randn()
    xyz[:,1] = xyz[:,1] + y + watdel2 + 0.2*random.randn()
    xyz[:,2] = xyz[:,2] + z1 + watdel
    new1.update_xyz(xyz)
    new1.serial = nres - 1
    new2 = pdb.Residue()
    new2.copy(watbox.residues[idx2]) #
    xyz = new2.collect("xyz")
    xyz = xyz - xyz[0,:]
    xyz[:,0] = xyz[:,0] + x + watdel2 + 0.2*random.randn()
    xyz[:,1] = xyz[:,1] + y + watdel2 + 0.2*random.randn()
    xyz[:,2] = xyz[:,2] + z2 - watdel
    new2.update_xyz(xyz)
    new2.serial = nres - 1
    membrane.residues.append(new1)
    membrane.residues.append(new2)
    y = y + watdel
    if y > memxyz[:,1].max() :
      y = memxyz[:,1].min()
      x = x + watdel
      if x > memxyz[:,0].max() :
        x = memxyz[:,0].min()
        z1 = z1 + watdel
        z2 = z2 - watdel

  natom = 0
  memxyz = np.zeros([nres*len(lipid.residues[0].atoms)+toadd*(len(watbox.residues[0].atoms)),3])
  for residue in membrane.residues :
    for atom in residue.atoms :
      natom = natom + 1
      atom.serial = natom
      atom.residue = residue.serial
      memxyz[natom-1,:] = atom.xyz
  membrane.write(args.out)
  len = memxyz.max(axis=0)-memxyz.min(axis=0)
  print "X = %.3f Y = %.3f Z = %.3f"%(len[0],len[1],len[2])
  print "editconf -f %s -o %s -box %.3f %.3f %.3f"%(args.out,args.out[:-3]+"gro",len[0]/10.0,len[1]/10.0,len[2]/10.0)
  
  
  
  
  
