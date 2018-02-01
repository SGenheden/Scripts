# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to perform actions on MD trajectories
"""

import os
from collections import namedtuple

import MDAnalysis as md
import MDAnalysis.core.AtomGroup as AtomGroup
import MDAnalysis.analysis.align as align
import MDAnalysis.lib.util as mdutil
import numpy as np
try :
    import pyvoro
except:
    pass
from scipy.spatial.distance import cdist
import sklearn.mixture as mixture

from . import pbc
from . import geo
from . import groups
from . import mol

MDRecord = namedtuple("MDRecord",["time","value"])

resanallib = {
    "me" : {
        'VAL' : [('CB','CG1'),('CB','CG2')],
        'THR' : [('CB','CG2')],
        'ILE' : [('CB','CG2'),('CG1','CD1')],
        'LEU' : [('CG','CD1'),('CG','CD2')],
        'MET' : [('SD','CE')],
        'ALA' : [('CA','CB')]
    },
    "dict" : {
        'VAL' : [("CB","CG1"),("CB","CG2")],
        'SER' : [("CB","HB2"),("CB","HB3")],
        'THR' : [("CB","CG2")],
        'ILE' : [("CG1","CD1")],
        'LEU' : [("CG","CD1"),("CG","CD2")],
        'MET' : [("SD","CE")],
        'ASN' : [("ND2","HD21"),("ND2","HD22")],
        'GLN' : [("NE2","HE21"),("NE2","HE22")],
        'PHE' : [("CD1","HD1")],
        'HID' : [("CD2","HD2")],
        'HIE' : [("CD2","HD2")],
        'HIP' : [("CD2","HD2")],
        'TYR' : [ ("CD1","HD1")],
        'PRO' : [("CG","HG2"),("CG","HG3")],
        'LYS' : [("CB","HB2"),("CB","HB3"),("CG","HG2"),("CG","HG3"),("CD","HD2"),("CD","HD3"),("CE","HE2"),("CE","HE3")],
        'ARG' : [("CB","HB2"),("CB","HB3"),("CG","HG2"),("CG","HG3"),("CD","HD2"),("CD","HD3"),("NE","HE")],
        'ASP' : [("CB","HB2"),("CB","HB3")],
        'GLU' : [("CG","HG2"),("CG","HG3")],
    }
}

class TrajectoryAction(object):
    """
    Base class for actions that are called by a TrajectoryProcessor
    object. Classes that wants to implement a particular action
    should inherit from this.
    """
    def __init__(self, processor):
        self.processor = processor
        self.dosubsample = False
        processor.actions.append(self)

    def add_arguments(self):
        """
        Function that is called to add action-specific arguments

        Arguments
        ---------
        parser : argparse.ArgumentParser
            the parser to add command-line arguments to
        """
        pass

    def setup(self, args):
        """
        Function that is called after the processor has parsed the command-line
        arguments.

        Arguments
        ---------
        args : argparse.Namespace
            the parsed arguments
        """

    def process(self):
        """
        Main processing function called at each trajectory timestep
        """
        pass

    def subsample(self):
        """
        Function called occasionally if subsampling is turned on
        """
        pass

    def finalize(self):
        """
        Function that is called after all timesteps has been processed
        """
        pass

    def _write_records(self, postfix='', headers=None):
        """
        Helper routine to write out a list of MDRecord to disc
        """
        if self.records :
            with open(self.out+postfix,'w') as f :
                if headers is not None :
                    f.write("#"+"\t".join(headers)+"\n")
                for entry in self.records:
                    if isinstance(entry.value,float) or isinstance(entry.value,np.float32):
                        f.write("%.0f\t%.3f\n"%(entry.time,entry.value))
                    elif isinstance(entry.value,int) or isinstance(entry.value,np.int32):
                        f.write("%.0f\t%d\n"%(entry.time,entry.value))
                    else:
                        f.write("%.0f\t%s\n"%(entry.time,"\t".join("%.3f"%v for v in entry.value)))


ResidueAtoms = namedtuple("ResidueAtoms",["first","last"])

class AnalysisGrid(object) :

    """
    Class to store a 2D grid for analysis

    Attributes
    ----------
    count : NdArray
        the number of times the accumulate function has been called
    edgesx : NdArray
        the edges along the x-dimension
    edgesy : NdArray
        the edges along the y-dimension
    matrix : NdArray
        the discretised data
    resolution : float
        the resolution, i.e. the size of each pixel
    """
    def __init__(self, xyz, resolution=1.0) :
        xyz_shift = xyz - xyz.mean(axis=0)

        self.edgesx = np.arange(xyz_shift[:,0].min(),
                                    xyz_shift[:,0].max()+resolution, resolution)
        self.edgesy = np.arange(xyz_shift[:,1].min(),
                                    xyz_shift[:,1].max()+resolution, resolution)
        self.resolution = resolution
        self.matrix = np.zeros([self.edgesx.shape[0]+1, self.edgesy.shape[0]+1])
        self.count = np.zeros([self.edgesx.shape[0]+1, self.edgesy.shape[0]+1])

    def accumulate(self, positions, data) :

        idx = self.indices(positions)
        self.matrix[idx[0], idx[1]] += data
        self.count[idx[0], idx[1]] += 1.0

    def average(self, func=None) :

        sel = self.count > 0.0
        self.matrix[sel] /= self.count[sel]
        if func is not None :
            self.matrix[sel] = func(self.matrix[sel])

    def indices(self, positions) :
        """
        Returns the grid coordinates for a set of Cartesian coordinates
        """
        xidx = np.digitize(positions[:,0], self.edgesx)
        yidx = np.digitize(positions[:,1], self.edgesy)
        return xidx, yidx

    def write(self, filename) :

        with open(filename, "w") as f :
            for i in range(self.matrix.shape[0]) :
                for j in range(self.matrix.shape[1]) :
                    f.write("%.3f\t"%self.matrix[i,j])
                f.write("\n")

class CenterWholeAlign(TrajectoryAction):
    """
    Class to make MD snapshots whole over periodic boxes and to centre and
    align proteins.

    Attributes
    ----------
    protsel : MDAnalysis.AtomGroup
        the protein selection
    refuni : MDAnalyis.Universe
        the reference universe used for alignment
    residue_atoms : list of integer tuples
        the atom numbers of all residues excluding proteins
    residues : list of MDAnalysis.AtomGroup
        the residues in the universe excluding proteins
    records : list of MDRecords
        the RMSD record at each processes snapshot
    writer : MDAnalysis.Writer
        the output trajectory writer
    """
    def add_arguments(self, parser):
        parser.add_argument('--bbmask',help="the selectiom mask for backbone",default="name CA")
        parser.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
        parser.add_argument('-o','--out',help="the output",default="centerwhole")
        parser.add_argument('--noalign',action="store_true",help="turns off alignment",default=False)
        parser.add_argument('--nocenter',action="store_true",help="turns off centering",default=False)
        parser.add_argument('--nowhole',action="store_true",help="turns off making whole",default=False)

    def setup(self, args):
        self.refuni = md.Universe(self.processor.args.struct)
        self.protsel = self.processor.universe.select_atoms(args.pmask)
        if len(self.protsel) == 0 :
            self.nocenter = True
            self.noalign = True
        else:
            self.nocenter = args.nocenter
            self.noalign = args.noalign
        self.nowhole = args.nowhole

        self.residues = []
        self.residue_atoms = []
        for res in self.processor.universe.select_atoms("not "+args.pmask).residues:
            if len(res.atoms) > 1 :
                self.residues.append(res)
            self.residue_atoms.append(ResidueAtoms(res.atoms[0].index,res.atoms[-1].index))

        self.records = []
        self.writer = md.Writer(args.out,
                                    self.processor.universe.trajectory.n_atoms)
        self.out = args.out
        self.bbmask = args.bbmask

    def process(self):
        if not self.nowhole :
            if len(self.protsel) > 0:
                xyz = pbc.make_whole_xyz(self.protsel.positions,self.processor.currbox)
                self.protsel.positions = xyz
            for res in self.residues :
                xyz = pbc.make_whole_xyz(res.atoms.positions,self.processor.currbox)
                res.atoms.positions = xyz
        if not self.nocenter :
            self._center()
        if not self.noalign :
            rmsd = align.alignto(self.processor.universe, self.refuni,
                                select=self.bbmask)[1]
            self.records.append(MDRecord(self.processor.currtime,rmsd))

        self.writer.write(self.processor.currsnap)

    def finalize(self):
        """
        Write out the RMSDs to disc and close the output trajectory
        """
        self._write_records(postfix="_rmsd.txt")
        try :
            self.writer.close_trajectory()
        except :
            pass

    def _center(self) :
        xyz = self.processor.currsnap._pos
        #com1 = xyz[self.protsel[0].number:self.protsel[-1].number+1].mean(axis=0)
        com1 = self.protsel.center_of_geometry()

        for residue in self.residue_atoms :
            com2 = xyz[residue.first:residue.last+1].mean(axis=0)
            dr = pbc.unwrap_vector(com1 - com2, self.processor.currbox)
            xyz[residue.first:residue.last+1] = xyz[residue.first:residue.last+1] + dr

        delta = com1 - self.processor.currbox/2.0
        self.processor.currsnap._pos = xyz - delta

MDGroupSelection = namedtuple("MDGroupSelection",["atomgroup", "indices", "transmat"])

class ChainOrderAnalysis(TrajectoryAction):
    """
    Class to analyse chain order parameters during a trajectory

    Attributes:
    -----------
    normal : numpy ndarray
      the normal of the membrane, assumed to be z-axis
    out : string
      the output filename
    selections : list of MDAnalysis.AtomGroup
        the  selections
    records : list of MDRecord
        the chain orders at each timestep
    """
    def add_arguments(self, parser):
        parser.add_argument('--selections',nargs="+", help="the chains")
        parser.add_argument('--analysis',choices=["CC","CH"], help="the type of analysis C-C or C-H", default="CC")
        parser.add_argument('--groups', help="group definitions for pseudo-atom calculation")
        parser.add_argument('--gridout', help="the prefix for the filename of a 2D grid")
        parser.add_argument('--protmask',help="the selectiom mask for lipid residues")
        parser.add_argument('--pmask',help="the selectiom mask for phosphor atoms",default="name P")
        parser.add_argument('-o', '--out', help="the output prefix", default="order")

    def setup(self, args):

        def _get_h(atomgrp):
            for atom2 in atomgrp[0].bonded_atoms :
                if atom2.mass < 5.0 :
                    return self.processor.universe.select_atoms("resname %s and name %s"%(atomgrp[0].resname,atom2.name))
            raise Exception("Could not find any H atom bonded to %s in %s"%(atomgrp[0].name,atomgrp[0].resname))

        def _enumerateatoms(resname, atomstr) :

            lipid = self.processor.universe.select_atoms("resname %s"%resname)[0].residue
            base = atomstr[:-1]
            atomi = int(atomstr[-1])
            lst = []
            while True :
                try :
                    name = "%s%d"%(base,atomi)
                    dummy = lipid[name]
                    lst.append(name)
                    atomi += 1
                except :
                    break
            return lst

        def _expandlist(liststr):

            l, r = liststr.split("..")

            i = l.find("(")
            start = int(l[i+1:])
            l = l[:i]

            i = r.find(")")
            end = int(r[:i])
            r = r[i+1:]

            return ["%s%d%s"%(l,i,r) for i in range(start,end+1)]

        self.headers = ["Time"]
        self.selheaders = []
        self.selections = []
        self.analtype = args.analysis
        if self.analtype == "CH" :
            self.hselections = []

        self.resgroups = None
        if args.groups is not None:
            if self.analtype == "CH" :
                raise Exception("Cannot perform C-H analysis on pseudo-atoms")
            self.resgroups = groups.read_groups(args.groups)

        for selin in args.selections:
            resname, chainlist = selin.split(":")

            if self.resgroups is not None:
                if resname not in self.resgroups:
                    raise Exception("Cannot find %s in groups spec."%resname)
                pseudoatoms = [group.name for group in self.resgroups[resname].groups]

            if chainlist.find("-") > -1:
                atomlist = chainlist.split("-")
            elif chainlist.find("..") > -1:
                atomlist = _expandlist(chainlist)
            elif chainlist.startswith("@") :
                atomlist = _enumerateatoms(resname, chainlist[1:])
            else:
                raise Exception("Atom list need be specified with '-' or with expansion '(..)'")

            if self.resgroups is None:
                atomsels = [self.processor.universe.select_atoms("resname %s and name %s"%(resname,atom))
                        for atom in atomlist]
                print "%s (%s) - %d atoms and %d atoms in first selection"% \
                    (resname, ",".join(atomlist), len(atomlist), len(atomsels[0]))
                for atomgrp, atom in zip(atomsels[1:], atomlist[1:]):
                    if len(atomgrp) != len(atomsels[0]):
                        raise Exception("Selection for %s is different in length than the first selection"%atom)
                self.selections.append(atomsels)
            else:
                for atom in atomlist:
                    if atom not in pseudoatoms :
                        raise Exception("Could not find selected atom %s in the group spec."%atom)
                # Select all atoms for the selected residue, the coordinates
                # will be transformed to pseudo-atoms
                atomsel = self.processor.universe.select_atoms("resname %s"%resname)
                atomnames = [atom.name for atom in atomsel.residues[0].atoms]
                ngroups = len(self.resgroups[resname].groups)
                natoms = len(atomnames)
                nres = len(atomsel.residues)
                # Create the pseudo atom indices
                indices0 = [self.resgroups[resname].indices(atom) for atom in atomlist]
                indices = [[i0[0]+ngroups*i for i in range(nres)] for i0 in indices0]
                # Create the transformation matrix by replacting the one for the first residue
                transmat0 = self.resgroups[resname].transmat(atomnames)
                transmat = np.zeros([ngroups*nres,natoms*nres])
                for i in range(nres):
                    transmat[i*ngroups:(i+1)*ngroups,i*natoms:(i+1)*natoms] = transmat0
                self.selections.append(MDGroupSelection(atomsel, indices, transmat))
                print "%s (%s) - %d atoms and %d atoms in first selection"% \
                    (resname, ",".join(atomlist), len(atomlist), len(indices[0]))

            self.headers.extend(["%s/%s"%(resname, atom) for atom in atomlist])
            self.selheaders.append(["%s/%s"%(resname, atom) for atom in atomlist])

            if self.analtype == "CH":
                hatomsels = [_get_h(atomgrp) for atomgrp in atomsels]
                self.hselections.append(hatomsels)
                for atomgrp, atom in zip(hatomsels[1:], atomlist[1:]):
                    if len(atomgrp) != len(hatomsels[0]):
                        raise Exception("H-selection for %s is different in length than the first selection"%atom)

        self.out = args.out
        # Assumes that the normal is along the z-axis
        self.normal = np.array([0.0,0.0,1.0])
        self.records = []

        self.gridout = args.gridout
        self.phosphorsel = None
        if self.gridout is not None :
            bounds = np.asarray([[0.0, 0.0, 0.0],self.processor.universe.dimensions[:3]])
            self.grid_low = AnalysisGrid(bounds)
            self.grid_upp = AnalysisGrid(bounds)
            self.phosphorsel = self.processor.universe.select_atoms(args.pmask)
            if args.protmask is not None :
                self.protsel = self.processor.universe.select_atoms(args.protmask)
                self.grid_prot = AnalysisGrid(bounds)
                self.protone = np.ones(len(self.protsel))
            else :
                self.protsel = None

    def process(self):

        mid = None
        if self.gridout is not None :
            mid = self.phosphorsel.center_of_geometry()
            if self.protsel is not None :
                self.grid_prot.accumulate(self.protsel.positions-mid, self.protone)

        orders = []
        if self.analtype == "CC":
            if self.resgroups is None:
                for selection in self.selections :
                    for a1, a2 in zip(selection[:-1],selection[1:]):
                        orders.append(self._calc_order(a1.positions,
                                        a2.positions, self.normal, mid))
            else:
                if self.processor.nprocessed == 1:
                    f = open(self.out+"_first_pseudo.xyz", "w")
                for selection in self.selections :
                    xyz = np.dot(selection.transmat, selection.atomgroup.positions)
                    if self.processor.nprocessed == 1:
                        for pos in xyz:
                            f.write("c %.3f %.3f %.3f\n"%(pos[0], pos[1], pos[2]))
                    for i1, i2 in zip(selection.indices[:-1], selection.indices[1:]):
                        orders.append(self._calc_order(xyz[i1,:],
                                        xyz[i2,:], self.normal, mid))
                if self.processor.nprocessed == 1:
                    f.close()
        elif self.analtype == "CH":
            for cselection, hselection in zip(self.selections, self.hselections):
                for a1, a2 in zip(cselection, hselection):
                    orders.append(self._calc_order(a1.positions,
                                    a2.positions, self.normal, mid))
        self.records.append(MDRecord(self.processor.currtime, orders))

    def _calc_order(self, a1, a2, norm, mid):
        # Atom2 - Atom1
        vec = a2 - a1
        # Projection with normal
        proj = np.multiply(vec,norm).sum(axis=1)**2 / np.sum(vec**2,axis=1)

        # Discretize on a grid
        if self.gridout is not None :
            sel_low = self.phosphorsel.positions[:,2] < mid[2]
            sel_upp = np.logical_not(sel_low)
            coords_upp = self.phosphorsel.positions[sel_upp,:]
            coords_low = self.phosphorsel.positions[sel_low,:]
            self.grid_low.accumulate(coords_low-mid, proj[sel_low])
            self.grid_upp.accumulate(coords_upp-mid, proj[sel_upp])

        # return order parameter
        return np.abs(0.5*(3.0*proj.mean()-1))

    def finalize(self):
        self._write_records(postfix="_dt.txt", headers=self.headers)

        data = np.asarray([r.value for r in self.records])
        av = data.mean(axis=0)
        std = data.std(axis=0)
        offset = 0
        selavs = []
        selstds = []
        fac = -1 if self.analtype == "CC"  else 0
        for heads in self.selheaders:
            selavs.append(av[offset:offset+len(heads)+fac])
            selstds.append(std[offset:offset+len(heads)+fac])
            offset += len(heads)+fac
        maxatm = max([len(heads) for heads in self.selheaders])+fac

        with open(self.out+".txt", "w") as f :
            f.write("".join(["\t%s\t\t"%heads[0].split("/")[0] for heads in self.selheaders])+"\n")
            for i in range(maxatm):
                for j in range(len(self.selheaders)):
                    if i < len(self.selheaders[j]) :
                        f.write("%s\t%.3f\t%.3f\t"%(self.selheaders[j][i].split("/")[1],
                            selavs[j][i],selstds[j][i]))
                    else:
                        f.write(" \t \t \t")
                f.write("\n")
            for avs in selavs:
                f.write("Av\t%.3f\t%.3f\t"%(avs.mean(),avs.std()/np.sqrt(avs.shape[0])))
            f.write("\n")

        if  self.gridout is not None :
            def order(mat) :
                return np.abs(0.5*(3.0*mat-1))
            self.grid_low.average(func=order)
            self.grid_low.write(self.gridout+"_low.dat")
            self.grid_upp.average(func=order)
            self.grid_upp.write(self.gridout+"_upp.dat")
            if self.protsel is not None :
                self.grid_prot.average()
                self.grid_prot.write(self.gridout+"_prot.dat")

class IredAnalysis(TrajectoryAction):
    """
    Analysis class for iRED

    Attributes
    ----------
    atm1 : MDAnalysis.AtomSelection
        the first atom making up the iRED vector
    atm2 : MDAnalysis.AtomSelection
        the second atom making up the iRED vector
    mat : numpy.ndarray
        the built-up of the correlation matrix
    s2list : list
        the S2 order parameter for each vector at each subsample point
    outname : string
        the name of the output file
    processor : TrajectoryProcessor object
        the trajectory processor calling this analysis
    """
    def add_arguments(self, parser):
        parser.add_argument('--atoms',nargs=2,help="the atom names making the vectors",default=["N","H"])
        parser.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
        parser.add_argument('-o','--out',help="the output name",default="s2.txt")
        parser.add_argument('--resoffset',type=int,help="the residue offset",default=0)
        parser.add_argument('--uselib',choices=["no","dict","me"],help="if to use library vectors",default="no")


    def setup(self,args):
        protsel = self.processor.universe.select_atoms(args.pmask)
        self.uselib = args.uselib
        if self.uselib == "no":
            self.atm2 = protsel.select_atoms("name "+args.atoms[1])
            self.atm1 = protsel.select_atoms("name "+args.atoms[0]+
                                        " and byres name "+args.atoms[1])
        elif self.uselib in resanallib :
            lib = resanallib[self.uselib]
            atm1 = []
            atm2 = []
            for res in self.processor.universe.select_atoms("protein").residues:
                if res.name not in lib :
                    continue
                for atompairs in lib[res.name]:
                    atm1.append(res[atompairs[0]])
                    atm2.append(res[atompairs[1]])
            self.atm2 = AtomGroup.AtomGroup(atm2)
            self.atm1 = AtomGroup.AtomGroup(atm1)
        self.mat = np.zeros([len(self.atm1),len(self.atm1)])
        self.s2list = []
        self.outname = args.out
        self.dosubsample = True
        self.resoffset = args.resoffset

    def process(self):
        """
        Building up the correlation matrix, called at each MD snapshot
        """
        v1 = self.atm2.positions-self.atm1.positions
        vlen = 1.0 / np.sqrt((v1*v1).sum(axis=1))

        """mat2 = np.zeros(mat.shape)
        for i in range(nvec):
            for j in range(nvec):
                mat2[j,i] = np.sum(v1[i]*v1[j])*(vlen[i]*vlen[j]) """

        xx1,xx2 = np.meshgrid(v1[:,0],v1[:,0])
        yy1,yy2 = np.meshgrid(v1[:,1],v1[:,1])
        zz1,zz2 = np.meshgrid(v1[:,2],v1[:,2])
        ll1,ll2 = np.meshgrid(vlen,vlen)
        mat0 = (xx1*xx2+yy1*yy2+zz1*zz2)*(ll1*ll2)
        self.mat += 3.0*mat0*mat0-1

    def subsample(self):
        """
        Calculating the S2 order parameters and then zero the correlation matrix
        """
        self.mat = 0.5*(self.mat / float(self.processor.subsamples))
        # Calculating and sorting the eigenvalues and eigenvectors
        eval,evec = np.linalg.eig(self.mat)
        idx = eval.argsort()[::-1]
        eval = eval[idx]
        evec = evec[:,idx]

        prod = evec*evec
        s2 = np.array([1.0-np.sum(eval[5:]*prod[i,5:]) for i in range(prod.shape[0])])

        self.s2list.append(s2)
        self.mat = np.zeros(self.mat.shape)

    def finalize(self):
        """
        Write out the order parameters to disc
        """
        with open(self.outname,"w") as f :
            self.s2list = np.asarray(self.s2list)
            frmstr = "%.5f "*self.s2list.shape[0] # ///frmstr%tuple(rs2)
            prevres = None
            prevatom = None
            reslist = []
            for i,(atm,rs2) in enumerate(zip(self.atm1,self.s2list.T)):
                if self.uselib != "dict":
                    f.write("%s %d %.5f\n"%(atm.resname,atm.resnum+self.resoffset,rs2.mean()))
                    #f.write("%s %d %s\n"%(atm.resname,atm.resnum+self.resoffset," ".join("%.5f"%v for v in rs2)))
                else :
                    if reslist and atm.resnum != prevres:
                        av = np.asarray(reslist).mean()
                        f.write("%s %d %.5f\n"%(prevatom.resname,prevatom.resnum+self.resoffset,av))
                        reslist = []
                    reslist.append(rs2.mean())
                    prevres = atm.resnum
                    prevatom = atm
            if self.uselib == "dict":
                av = np.asarray(reslist).mean()
                f.write("%s %d %.5f\n"%(atm.resname,atm.resnum+self.resoffset,av))

class MemBulkAnalysis(TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--pmask',help="the selectiom mask for phosphor atoms",default="name P")
        parser.add_argument('--wmask',help="the selectiom mask for water atoms",default="name OH2")
        parser.add_argument('--smask',help="the selectiom mask for solute atoms")
        parser.add_argument('--sconc',type=float, help="the target solute concentration",default=1.0)
        parser.add_argument('--svol',type=float, help="the solute number volume",default=1.0)
        parser.add_argument('--wvol',type=float, help="the water number volume",default=0.0181)

    def setup(self, args):
        self.phosphorsel = self.processor.universe.select_atoms(args.pmask)
        self.watersel = self.processor.universe.select_atoms(args.wmask)
        self.allsel = self.processor.universe.select_atoms("all")
        print "Number of phosphor (%d) and water (%d) atoms"%(
                len(self.phosphorsel), len(self.watersel))

        if args.smask is not None :
            self.solute = self.processor.universe.select_atoms(args.smask)
            print "Number of solute atoms = %d"%len(self.solute)
        else :
            self.solute = None

        self.nphosphor = 1.0 / float(len(self.phosphorsel))
        self.nwater = 1.0 / float(len(self.watersel))

        # Setup edges to cover the entire simulation box
        zpos = self.processor.universe.coord._pos[:,2] - self.allsel.positions[:,2].mean()
        self.resolution = 0.25
        self.edges = np.arange(zpos.min(),zpos.max()+self.resolution,self.resolution)
        self.zvals = 0.5 * (self.edges[:-1] + self.edges[1:]) * 0.1

        self.pdensity = np.zeros(self.edges.shape[0]-1)
        self.wdensity = np.zeros(self.edges.shape[0]-1)
        self.wdensity_now = np.zeros(self.edges.shape[0]-1)
        if self.solute is not None :
            self.sdensity_now = np.zeros(self.edges.shape[0]-1)

        self.sconc = args.sconc
        self.svol = args.svol
        self.wvol = args.wvol

    def process(self) :

        zpos = self.phosphorsel.positions[:,2] - self.allsel.positions[:,2].mean()
        hist, b = np.histogram(zpos, bins=self.edges)
        self.pdensity += hist

        zpos = self.watersel.positions[:,2] - self.allsel.positions[:,2].mean()
        hist, b = np.histogram(zpos, bins=self.edges)
        self.wdensity += hist
        self.wdensity_now = hist

        if self.solute is not None :
            zpos = self.solute.positions[:,2] - self.allsel.positions[:,2].mean()
            hist, b = np.histogram(zpos, bins=self.edges)
            self.sdensity_curr = hist

    def finalize(self):

        # Calculate how many water molecules in the bulk
        firsti, lasti = mol.density_intercept(self.wdensity, self.pdensity)
        nbulkwat = int(np.round(self.wdensity_now[:firsti].sum()
                                +self.wdensity_now[lasti+1:].sum()))
        print "Nwat\t%d"%nbulkwat

        # Calculate how many solutes there are outside the membrane
        if self.solute is not None :
            nbulksol = int(np.round(self.sdensity_now[:firsti].sum()
                                    +self.sdensity_now[lasti+1:].sum()))
        else :
            nbulksol = 0
        print "Nsol\t%d"%nbulksol

        ntot = nbulkwat + nbulksol
        ntarget = np.round(ntot*self.wvol/(1.0/self.sconc+self.wvol-self.svol))
        nadd = ntarget-nbulksol
        if nadd % 2 == 0 :
            nadd = (0.5*nadd , 0.5*nadd)
        else :
            nadd = (np.ceil(0.5*nadd) , np.floor(0.5*nadd))
        print "Nadd\t%d\t%d"%nadd



class MemDensAnalysis(TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--pmask',help="the selectiom mask for phosphor atoms",default="name P")
        parser.add_argument('--wmask',help="the selectiom mask for water atoms",default="name OH2")
        parser.add_argument('--smask',help="the selectiom mask for solute atoms")
        parser.add_argument('-o','--out',help="the output prefix",default="memdens")

    def setup(self, args):
        self.dosubsample = True
        self.out = args.out
        self.phosphorsel = self.processor.universe.select_atoms(args.pmask)
        self.watersel = self.processor.universe.select_atoms(args.wmask)
        self.allsel = self.processor.universe.select_atoms("all")
        print "Number of phosphor (%d) and water (%d) atoms"%(
                len(self.phosphorsel), len(self.watersel))

        if args.smask is not None :
            self.solute = self.processor.universe.select_atoms(args.smask)
            print "Number of solute atoms = %d"%len(self.solute)
        else :
            self.solute = None

        self.nphosphor = 1.0 / float(len(self.phosphorsel))
        self.nwater = 1.0 / float(len(self.watersel))

        # Setup edges to cover the entire simulation box
        zpos = self.processor.universe.coord._pos[:,2] - self.allsel.positions[:,2].mean()
        self.resolution = 0.25
        self.edges = np.arange(zpos.min(),zpos.max()+self.resolution,self.resolution)
        self.zvals = 0.5 * (self.edges[:-1] + self.edges[1:]) * 0.1

        self.ppos_curr = []
        self.pdensity_curr = np.zeros(self.edges.shape[0]-1)
        self.wdensity_curr = np.zeros(self.edges.shape[0]-1)
        self.pdensity = []
        self.wdensity = []
        if self.solute is not None :
            self.sdensity_curr = np.zeros(self.edges.shape[0]-1)
            self.solute_snapshots = 0
            self.sdensity = []

        self.records = []

    def process(self) :

        zpos = self.phosphorsel.positions[:,2] - self.allsel.positions[:,2].mean()
        hist, b = np.histogram(zpos, bins=self.edges)
        self.pdensity_curr += hist
        self.ppos_curr.extend(zpos)

        zpos = self.watersel.positions[:,2] - self.allsel.positions[:,2].mean()
        hist, b = np.histogram(zpos, bins=self.edges)
        self.wdensity_curr += hist

        if self.solute is not None :
            zpos = self.solute.positions[:,2] - self.allsel.positions[:,2].mean()
            hist, b = np.histogram(zpos, bins=self.edges)
            self.sdensity_curr += hist
            self.solute_snapshots += 1

    def subsample(self) :

        # Calculate D_hh
        model = mixture.GaussianMixture(n_components=2)
        model.fit(np.asarray(self.ppos_curr).reshape(-1,1))
        dhh  = np.abs(model.means_[1][0]-model.means_[0][0]) * 0.1
        self.ppos_curr = []

        # Calculate intercept of water and phosphor density,
        # and from that the membrane volume
        firsti, lasti = mol.density_intercept(self.wdensity_curr, self.pdensity_curr)
        firstz = self.zvals[firsti]
        lastz = self.zvals[lasti]
        memfrac = (lastz - firstz) / (self.processor.currbox[2] * 0.1)

        # Calculate how many solutes there are inside the membrane
        if self.solute is not None :
            nfreq = 1 / float(self.processor.freq)
            solute_dens = self.sdensity_curr / float(self.solute_snapshots)
            ninside = int(np.round(solute_dens[firsti:lasti+1].sum()))
            self.records.append(MDRecord(self.processor.currtime, [dhh, lastz - firstz,
                            self.processor.currbox[2] * 0.1, memfrac, ninside]))

        else :
            self.records.append(MDRecord(self.processor.currtime, [dhh, lastz - firstz,
                self.processor.currbox[2] * 0.1, memfrac]))

        # Store away the accumulayed densities and zero them
        self.pdensity.append(self.pdensity_curr)
        self.wdensity.append(self.wdensity_curr)
        self.pdensity_curr = np.zeros(self.edges.shape[0]-1)
        self.wdensity_curr = np.zeros(self.edges.shape[0]-1)

        if self.solute is not None :
            self.sdensity.append(self.sdensity_curr)
            self.sdensity_curr = np.zeros(self.edges.shape[0]-1)
            self.solute_snapshots = 0

    def finalize(self):

        def _write_density(density, scaling, postfix) :
            density = np.asarray(density) * scaling
            with open(self.out+postfix, "w") as f :
                for z, av, err in zip(self.zvals, density.mean(axis=0),
                        density.std(axis=0)/np.sqrt(density.shape[0])) :
                    f.write("%.3f %.3f %.3f\n"%(z, av, err))

        _write_density(self.pdensity, self.nphosphor, "_pdens.dat")
        _write_density(self.wdensity, self.nwater, "_wdens.dat")
        if self.solute is not None :
            _write_density(self.sdensity, 1.0 / len(self.solute), "_sdens.dat")
        self._write_records(postfix="_dt.txt")

        vals = np.asarray([entry.value for entry in self.records])
        with open(self.out+".txt", "w") as f :
            f.write(" ".join("%.3f %.3f"%(av, err) for av, err in zip(vals.mean(axis=0),
                        vals.std(axis=0)/np.sqrt(vals.shape[0])))+"\n")


class MempropAnalysis(TrajectoryAction):
    def add_arguments(self, parser):
        parser.add_argument('--pmask',help="the selectiom mask for phosphor atoms",default="name P")
        parser.add_argument('--lipidmask',help="the selectiom mask for lipid residues",default="resname POPC")
        parser.add_argument('--watmask',help="the selectiom mask for water residues",default="resname SOL")
        parser.add_argument('--watvol',type=float,help="the volume of a water molecule in nm3",default=0.0306)
        parser.add_argument('--gridout', help="the prefix for the filename of a 2D grid")
        parser.add_argument('--protmask',help="the selectiom mask for lipid residues")
        parser.add_argument('-o','--out',help="the output prefix",default="memprop")

    def setup(self,args):
        self.out = args.out
        self.phosphorsel = self.processor.universe.select_atoms(args.pmask)
        self.lipidsel = self.processor.universe.select_atoms(args.lipidmask)
        watsel  = self.processor.universe.select_atoms(args.watmask)
        self.nlipid = len(self.lipidsel.residues)
        self.nwat = len(watsel.residues)
        nphosph = len(self.phosphorsel.residues)
        print "Number of lipids (%d), waters (%d) and phosphor atoms (%d)"%(self.nlipid,self.nwat,nphosph)
        self.watvol = args.watvol
        if self.nlipid == 0 or self.nwat == 0 or nphosph == 0 :
            raise Exception("Either number of lipids (%d), water (%d) or phosphor atoms (%d)  is zero"%(self.nlipid,self.nwat,nphosph))
        self.apllist = []
        self.vpllist = []
        # Setup edges to cover the entire simulation box
        zpos = self.processor.universe.coord._pos[:,2] - self.lipidsel.positions[:,2].mean()
        self.resolution = 0.25
        self.edges = np.arange(zpos.min(),zpos.max()+self.resolution,self.resolution)
        self.density = np.zeros(self.edges.shape[0]+1)
        # Setup arrays for RMSF calculations
        self.sumcoords2 = np.zeros([nphosph,2])
        self.sumcoords = np.zeros([nphosph,2])
        self.records = []

        self.gridout = args.gridout
        if self.gridout is not None :
            bounds = np.asarray([[0.0, 0.0, 0.0],self.processor.universe.dimensions[:3]])
            self.grid_low = AnalysisGrid(bounds)
            self.grid_upp = AnalysisGrid(bounds)
            if args.protmask is not None :
                self.protsel = self.processor.universe.select_atoms(args.protmask)
                self.grid_prot = AnalysisGrid(bounds)
                self.protone = np.ones(len(self.protsel))
            else :
                self.protsel = None

    def process(self):
        """
        Calculate APL, VPL and accumulate density of phosphor selection
        """
        boxnm = self.processor.currbox / 10.0
        self.apllist.append(boxnm[0]*boxnm[1]/float(self.nlipid/2))
        self.vpllist.append((boxnm[0]*boxnm[1]*boxnm[2] -
                                self.watvol*self.nwat)/float(self.nlipid))
        zpos = self.phosphorsel.positions[:,2] - self.lipidsel.positions[:,2].mean()
        for lipdig in np.digitize(zpos,self.edges) :
            self.density[lipdig] += 1
        self.sumcoords += self.phosphorsel.positions[:,:2]
        self.sumcoords2 += self.phosphorsel.positions[:,:2]*self.phosphorsel.positions[:,:2]
        self.records.append(MDRecord(self.processor.currtime,[self.apllist[-1],self.vpllist[-1],self._calc_dhh(),self._calc_rmsf()]))

        if self.gridout is not None :
            mid = self.phosphorsel.center_of_geometry()
            sel_low = self.phosphorsel.positions[:,2] < mid[2]
            sel_upp = np.logical_not(sel_low)
            coords_upp = self.phosphorsel.positions[sel_upp,:]
            coords_low = self.phosphorsel.positions[sel_low,:]
            self.grid_low.accumulate(coords_low-mid,
                                        self._calc_zdist(coords_low, coords_upp))
            self.grid_upp.accumulate(coords_upp-mid,
                                        self._calc_zdist(coords_upp, coords_low))
            if self.protsel is not None :
                self.grid_prot.accumulate(self.protsel.positions-mid, self.protone)

    def finalize(self):
        """
        Calculate average APL and VPL as well as distance
        between peaks in the phosphor density
        """
        dhh = self._calc_dhh()
        apl = np.asarray(self.apllist).mean()
        vpl = np.asarray(self.vpllist).mean()
        rmsf = self._calc_rmsf()
        with open(self.out+".txt","w") as f :
            f.write("%.3f\t%.3f\t%.3f\t%.3f\n"%(apl, vpl, dhh, rmsf))
        self._write_records(postfix="_dt.txt")

        if  self.gridout is not None:
            self.grid_low.average()
            self.grid_low.write(self.gridout+"_low.dat")
            self.grid_upp.average()
            self.grid_upp.write(self.gridout+"_upp.dat")
            if self.protsel is not None :
                self.grid_prot.average()
                self.grid_prot.write(self.gridout+"_prot.dat")

    def _calc_dhh(self) :
        mid = int(self.density.shape[0]/2)
        dens_first = self.density[:mid]
        dens_last = self.density[mid:]
        max_first = np.argmax(dens_first)
        max_last = np.argmax(dens_last)
        return (max_last + mid - max_first) / 10.0 * self.resolution

    def _calc_rmsf(self):
        sumcoords = self.sumcoords / float(self.processor.nprocessed)
        sumcoords2 = self.sumcoords2 / float(self.processor.nprocessed)
        var = sumcoords2 - (sumcoords * sumcoords)
        return var.sum(axis=1).mean()*0.01

    def _calc_zdist(self, coords1, coords2) :
        """
        Calculate the z-distance between all lipids in one leaflet and the closest lipid in the other leaflet
        """
        dist = cdist(coords1[:,:2],coords2[:,:2],'sqeuclidean')
        j = np.argmin(dist,axis=1)
        return np.sqrt((coords2[j,2]-coords1[:,2])**2)*0.1

class MemVoronoiAnalysis(TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--mask',nargs="+",help="the selectiom mask for the atoms to do analysis on")
        parser.add_argument('--head',help="the name of the atom to determine leaflets",default="PO4")
        parser.add_argument('-o','--out',help="the output",default="memvoro")

    def setup(self, args):
        self.atoms = self.processor.universe.select_atoms(
            " or ".join("(%s)"%m for m in args.mask))
        self.head = self.processor.universe.select_atoms("name %s"%args.head)
        self.out = args.out
        self.aplrecords = []
        self.neighrecords = []
        self.resnames = list(set([atom.resname for atom in self.atoms]))
        self.respairs = []
        for i, resname1 in enumerate(self.resnames):
            for resname2 in self.resnames[i:]:
                self.respairs.append(resname1+"-"+resname2)

    def process(self):

        midz = self.head.positions[:,2].mean()
        lowsel = self.atoms.positions[:,2] < midz
        uppsel = np.logical_not(lowsel)

        celldim = [[0.0, self.processor.currbox[0]],
                    [0.0, self.processor.currbox[1]]]

        try :
            lareas, lneighbours = self._process_leaflet(self.atoms[lowsel], celldim)
            uareas, uneighbours = self._process_leaflet(self.atoms[uppsel], celldim)
        except:
            pass
        else:
            areas = 0.01 * 0.5 * (lareas + uareas)
            neighbours = 0.5 * (lneighbours + uneighbours)

            self.aplrecords.append(MDRecord(self.processor.currtime,areas))
            self.neighrecords.append(MDRecord(self.processor.currtime,neighbours))

    def _process_leaflet(self, atoms, celldim):

        cells = pyvoro.compute_2d_voronoi(atoms.positions[:,:2],
                celldim, 2.0, periodic=[True,True])


        # Calculate the area per each residue type
        areas = {resname : 0 for resname in self.resnames}
        nres  = {resname : 0.0 for resname in self.resnames}
        for atom, cell in zip(atoms, cells):
            areas[atom.resname] += cell["volume"]
            nres[atom.resname] += 1.0
        areaout = np.asarray([areas[resname] / nres[resname] for resname in self.resnames])

        # Calculate the neighbors
        vsets = [set((np.round(v[0],3),np.round(v[1])) for v in cell["vertices"]) for cell in cells]
        emptyset = set([])
        neighbors = {respair : 0 for respair in self.respairs}
        npairs = {respair : 0 for respair in self.respairs}
        for i, ivertices in enumerate(vsets):
            counts = {respair : 0 for respair in self.respairs}
            for j, jvertices in enumerate(vsets[i+1:],i+1):
                if ivertices & jvertices != emptyset :
                    iresname = atoms[i].resname
                    jresname = atoms[j].resname
                    if iresname+"-"+jresname in neighbors:
                        counts[iresname+"-"+jresname] += 1
                    else:
                        counts[jresname+"-"+iresname] += 1
            for respair in self.respairs:
                if counts[respair] > 0 :
                    npairs[respair] += 1.0
                    neighbors[respair] += counts[respair]
        neighout = np.asarray([neighbors[respair] / npairs[respair]
                                for respair in self.respairs])

        return areaout, neighout

    def finalize(self):
        headers = ["Time"]
        headers.extend(self.resnames)
        self.records = self.aplrecords
        self._write_records(postfix="_apl.txt", headers=headers)

        headers = ["Time"]
        headers.extend(self.respairs)
        self.records = self.neighrecords
        self._write_records(postfix="_neigh.txt", headers=headers)

class PrincipalAxisAnalysis(TrajectoryAction):
    """
    Class to analyse the principcal axis and its angle

    Attributes
    ----------
    masses : list of float
        the masses of the selected atoms
    normal : numpy.ndarray
        the normal to which the angle is calculate against
    records : list of MDRecord
        the recorded alpha (angle) values
    selection : MDAnalysis.AtomGroup
        the selection to make the analysis of
    """
    def add_arguments(self, parser):
        parser.add_argument('-m','--mask',help="the selectiom mask",default="name CA")
        parser.add_argument('-n','--normal',type=float,nargs=3,help="the normal vector",default=[0.0,0.0,1.0])
        parser.add_argument('-o','--out',help="the output filename",default="alpha.txt")

    def setup(self,args):
        self.selection = self.processor.universe.select_atoms(args.mask)
        self.masses = np.asarray([atom.mass for atom in self.selection])
        self.normal = np.asarray(args.normal)
        self.records = []
        self.out = args.out

    def process(self):
        xyz = pbc.make_whole_xyz(self.selection.positions,
                                    self.processor.currbox)
        moi = geo.moment_of_inertia(xyz-xyz.mean(axis=0),self.masses)
        princip = geo.principal_axes(moi)
        alpha = geo.angle(princip[0,:],self.normal)
        dalpha = pbc.unwrap_vector(alpha,np.pi)
        alpha = np.abs(alpha-dalpha)*180.0/np.pi
        self.records.append(MDRecord(self.processor.currtime,alpha))

    def finalize(self):
        """
        Write out average alpha and then all alphas to disc
        """
        alphas = np.asarray([entry.value for entry in self.records])
        print "Mean = %.3f Std = %.3f"%(alphas.mean(),alphas.std())
        self._write_records()

class RmsdAnalysis(TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--sel',help="the selectiom mask for the atoms to superpose",default="protein and (name C or name CA or name N or name O)")
        parser.add_argument('--mobile',help="the selectiom mask for the mobile atoms")
        parser.add_argument('-o','--out',help="the output",default="rmsd.out")

    def setup(self, args):
        self.refuni = md.Universe(self.processor.args.struct)
        self.sel = args.sel
        print "Will superpose %d atoms"%len(
            self.processor.universe.select_atoms(args.sel))
        self.out = args.out
        self.records = []
        if args.mobile is not None:
            self.refsel = self.refuni.select_atoms(args.mobile)
            self.trjsel = self.processor.universe.select_atoms(args.mobile)
            print "Mobile RMSD is on %d atoms"%len(self.refsel)
            self.mobrecords = []
            self.domobile = True
        else:
            self.domobile = False

    def process(self):
        rmsd = align.alignto(self.processor.universe, self.refuni,
                            select=self.sel)[1]
        self.records.append(MDRecord(self.processor.currtime,rmsd*0.1))

        if self.domobile :
            rmsd =  np.sqrt(np.mean((self.refsel.positions - self.trjsel.positions) ** 2))
            self.mobrecords.append(MDRecord(self.processor.currtime,rmsd*0.1))

    def finalize(self):
        self._write_records()
        if self.domobile:
            self.records = self.mobrecords
            r, t = os.path.splitext(self.out)
            self.out = r
            self._write_records(postfix="_mob"+t)

class RMSFAnalysis(TrajectoryAction):
    """
    Class to analyse the RMSF of selected residues

    Attributes
    ----------
    alignmask : string
        the mask to make the implicit RMSD fit on
    protsel : MDAnalysis.AtomGroup
        the selection used to select protein residues
    refuni : MDAnalysis.Universe
        the reference universe used for alignment
    sumcoords : numpy.ndarray
        the accumulated coordinates of the system
    sumcoords2 : numpy.ndarray
        the accumulated square of the system coordinates
    """
    def add_arguments(self, parser):
        parser.add_argument('--atoms',nargs="+",help="the atom names in the  backbone",default=["CA","N","C"])
        parser.add_argument('--fitsel',help="the selectiom mask for the atoms to superpose",default="protein and (name C or name CA or name N)")
        parser.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
        parser.add_argument('--uselib',choices=["no","dict","me"],help="if to use library vectors",default="no")
        parser.add_argument('--resoffset',type=int,help="the residue offset",default=0)
        parser.add_argument('-o','--out',help="the output",default="rmsf.txt")

    def setup(self, args):
        self.refuni = md.Universe(self.processor.args.struct)
        self.protsel = self.processor.universe.select_atoms(args.pmask)
        self.atoms = args.atoms
        self.fitsel = args.fitsel
        natm = len(self.processor.universe.atoms)
        self.sumcoords2 = np.zeros([natm,3])
        self.sumcoords = np.zeros([natm,3])
        self.out = args.out
        self.uselib = args.uselib
        self.resoffset = args.resoffset

    def process(self):
        rmsd = align.alignto(self.processor.universe, self.refuni, select=self.fitsel)
        xyz = self.processor.currsnap._pos
        self.sumcoords += xyz
        self.sumcoords2 += xyz*xyz

    def finalize(self):
        nsnap = float(self.processor.currtime / self.processor.dt)
        self.sumcoords = self.sumcoords / nsnap
        self.sumcoords2 = self.sumcoords2 / nsnap
        var = self.sumcoords2 - (self.sumcoords * self.sumcoords)
        bfac = (8.0/3.0)*np.pi*np.pi*var.sum(axis=1)*0.01 # Make the variance into nm^2

        with open(self.out,'w') as f:
            for residue in self.protsel.residues:
                if self.uselib == "no" :
                    rbfac = 0.0
                    masssum = 0.0
                    for atom in residue:
                        if atom.name not in self.atoms : continue
                        rbfac = bfac[atom.number]*atom.mass
                        masssum += atom.mass
                    f.write("%s %d %.3f\n"%(atom.resname,atom.resnum+self.resoffset,rbfac / masssum))
                elif self.uselib in resanallib:
                    lib = resanallib[self.uselib]
                    if residue.name not in lib:
                        continue
                    atomnames = [pair[0] for pair in lib[residue.name]]
                    for atom in residue:
                        if atom.name not in atomnames : continue
                        f.write("%s %d %.3f\n"%(atom.resname,atom.resnum+self.resoffset,bfac[atom.number]))

class SoluteGyration(TrajectoryAction):
    """
    Class to analyse radius of gyration for a solute

    Attributes
    ----------
    records : list of MDRecord
        the recorded alpha (angle) values
    selection : MDAnalysis.AtomGroup
        the selection to make the analysis of
    """
    def add_arguments(self, parser):
        parser.add_argument('-m','--mask',help="the selectiom mask")
        parser.add_argument('-o','--out',help="the output filename",default="gyration.txt")

    def setup(self,args):
        self.selection = self.processor.universe.select_atoms(args.mask).residues
        self.records = []
        self.out = args.out

    def process(self):
        g = [r.radius_of_gyration(pbc=False)*0.1 for r in self.selection]
        self.records.append(MDRecord(self.processor.currtime, g))

    def finalize(self):
        self._write_records()

class SoluteOrientation(TrajectoryAction):
    """
    Class to analyse the orientation of a solute

    Attributes
    ----------
    axis : ndarray
        the reference axis
    records : list of MDRecord
        the recorded alpha (angle) values
    sel : MDAnalysis.AtomGroup
        the selection mask
    """
    def add_arguments(self, parser):
        parser.add_argument('-m','--mask', nargs="+", help="the selectiom mask")
        parser.add_argument('--axis', nargs=3, type=float, help="the reference axis", default=[0.0,0.0,1.0])
        parser.add_argument('--zmask', help="the selection mask for z-partition")
        parser.add_argument('--zpart', nargs="+", type=float, help="the z-partition")
        parser.add_argument('-o','--out',help="the output filename",default="orient")

    def setup(self,args):
        if len(args.mask) == 1 :
            self.sel = self.processor.universe.select_atoms(args.mask).residues
            print "Selected %d"%len(self.sel)
        else :
            self.sel = [self.processor.universe.select_atoms(args.mask[0]),
                self.processor.universe.select_atoms(args.mask[1])]
            print "Selected %d"%len(self.sel[0])
        if args.zmask is not None :
            args.axis = [0.0,0.0,1.0]
            self.zmask = self.processor.universe.select_atoms(args.zmask)
            self.zpart = np.asarray(args.zpart)
        else:
            self.zmask = None
        self.axis = np.asarray(args.axis)
        self.records = []
        self.out = args.out

    def process(self):

        def _partz(i) :
            if self.zmask is None :
                return True
            else :
                zpos = self.zmask[i].position[2]
                if len(self.zpart) == 2 :
                    return zpos >= self.zpart[0] and zpos <= self.zpart[1]
                elif len(self.zpart) == 4 :
                    return (zpos >= self.zpart[0] and zpos <= self.zpart[1]) or \
                         (zpos >= self.zpart[2] and zpos <= self.zpart[3])
                else :
                    return True

        if len(self.sel) == 1 :
            vec = [r.principal_axes()[0] for i, r in enumerate(self.sel) if _partz(i)]
        else :
            vec = [a2.position - a1.position for i, (a1, a2) in enumerate(zip(self.sel[0],self.sel[1])) if _partz(i)]
        proj = np.asarray([np.dot(v, self.axis) / (np.linalg.norm(v)*np.linalg.norm(self.axis)) \
                            for v in vec])
        orient = np.arccos(proj)*180/np.pi
        self.records.append(MDRecord(self.processor.currtime, orient))

    def finalize(self):
        self._write_records(postfix="_dt.dat")

        if self.zmask is None :
            orient = np.asarray([r.value for r in self.records])
            avorient = orient.mean(axis=0)
            stdorient = orient.std(axis=0) #/ np.sqrt(orient.shape[0])
            proj = np.cos(orient)
            order = np.abs(0.5*(3.0*proj*proj-1).mean(axis=0))
            stdorder = 0.5*(3.0*proj*proj-1).std(axis=0) #/ np.sqrt(orient.shape[0])
            with open(self.out+".dat", "w") as f :
                f.write("\t".join(["%.3f\t%.3f"%(o,e) for o, e in zip(avorient, stdorient)]))
                f.write("\t")
                f.write("\t".join(["%.3f\t%.3f"%(o,e) for o, e in zip(order, stdorder)]))
                f.write("\n")
        with open(self.out+"_flat.dat", "w") as f :
            for r in self.records :
                for v in r.value :
                    f.write("%.3f\n"%v)

class StripAtoms(TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--sel',help="the selectiom mask for the atoms to strip")
        parser.add_argument('-o','--out',help="the output",default="stripped.dcd")

    def setup(self, args):
        self.sel = self.processor.universe.select_atoms("not ("+args.sel+")")
        print "Will save %d atoms to the new trajectory"%len(self.sel)
        notsel = self.processor.universe.select_atoms(args.sel)
        print "Will remove %d atoms"%len(notsel)
        self.writer = md.Writer(args.out,len(self.sel))
        self.out = args.out

    def process(self):
        self.writer.write(self.sel)

class VectorPlaneAngleAnalysis(TrajectoryAction):
    """
    Class to analyse the angle between a vector and a plane

    Attributes
    ----------
    normal : numpy.ndarray
        the normal to the plane
    records : list of MDRecord
        the recorded angle values
    atom1 : MDAnalysis.AtomGroup
        the selection of the first atom in the vector
    atom2 : MDAnalysis.AtomGroup
        the selection of the second atom in the vector
    """
    def add_arguments(self, parser):
        parser.add_argument('-a','--atoms',nargs=2,help="the selectiom mask for the two atoms")
        parser.add_argument('-n','--normal',type=float,nargs=3,help="the normal vector",default=[0.0,0.0,1.0])
        parser.add_argument('-o','--out',help="the output filename",default="angle.txt")

    def setup(self, args):
        self.atom1 = self.processor.universe.select_atoms(args.atoms[0])
        if len(self.atom1) == 0:
            raise Exception("Selection %s did not select any atoms"%args.atoms[0])
        self.atom2 = self.processor.universe.select_atoms(args.atoms[1])
        if len(self.atom2) == 0:
            raise Exception("Selection %s did not select any atoms"%args.atoms[1])
        elif len(self.atom2) != len(self.atom1):
            raise Exception("Selection %s and %s did not select the same number of atoms"%(args.atoms[0],args.atoms[1]))
        self.normal = np.asarray(args.normal)
        self.out = args.out
        self.records = []

    def process(self):
        vec = self.atom2.positions-self.atom1.positions
        sinang = np.multiply(vec,self.normal).sum(axis=1) / np.sqrt(np.sum(vec**2,axis=1))
        ang = np.rad2deg(np.arcsin(sinang))
        self.records.append(MDRecord(self.processor.currtime,ang))

    def finalize(self):
        allang = np.asarray([entry.value for entry in self.records])
        print "Mean = %.3f Std = %.3f"%(allang.mean(),allang.std())
        with open(self.out,'w') as f :
            for entry in self.records:
                f.write("%.0f %s\n"%(entry.time," ".join(["%.3f"%a for a in entry.value])))
