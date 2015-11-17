# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to perform actions on MD trajectories
"""

from collections import namedtuple

import MDAnalysis as md
import MDAnalysis.core.AtomGroup as AtomGroup
import MDAnalysis.analysis.align as align
import numpy as np

from . import pbc
from . import geo

MDRecord = namedtuple("MDRecord",["time","value"])

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

    def _write_records(self,postfix=''):
        """
        Helper routine to write out a list of MDRecord to disc
        """
        if self.records :
            with open(self.out+postfix,'w') as f :
                for entry in self.records:
                    if isinstance(entry.value,float) or isinstance(entry.value,np.float32):
                        f.write("%.0f\t%.3f\n"%(entry.time,entry.value))
                    elif isinstance(entry.value,int) or isinstance(entry.value,np.int32):
                        f.write("%.0f\t%d\n"%(entry.time,entry.value))
                    else:
                        f.write("%.0f\t%s\n"%(entry.time,"\t".join("%.3f"%v for v in entry.value)))


ResidueAtoms = namedtuple("ResidueAtoms",["first","last"])

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
        self.protsel = self.processor.universe.selectAtoms(args.pmask)
        if len(self.protsel) == 0 :
            self.nocenter = True
            self.noalign = True
        else:
            self.nocenter = args.nocenter
            self.noalign = args.noalign
        self.nowhole = args.nowhole

        self.residues = []
        self.residue_atoms = []
        for res in self.processor.universe.selectAtoms("not "+args.pmask).residues:
            if len(res) > 1 :
                self.residues.append(res)
            self.residue_atoms.append(ResidueAtoms(res[0].number,res[-1].number))

        self.records = []
        self.writer = md.Writer(args.out,
                                    self.processor.universe.trajectory.numatoms)
        self.out = args.out
        self.bbmask = args.bbmask

    def process(self):
        if not self.nowhole :
            if len(self.protsel) > 0:
                xyz = pbc.make_whole_xyz(self.protsel.get_positions(),self.processor.currbox)
                self.protsel.set_positions(xyz)
            for res in self.residues :
                xyz = pbc.make_whole_xyz(res.get_positions(),self.processor.currbox)
                res.set_positions(xyz)
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
        self.writer.close_trajectory()

    def _center(self) :
        xyz = self.processor.currsnap._pos
        #com1 = xyz[self.protsel[0].number:self.protsel[-1].number+1].mean(axis=0)
        com1 = self.protsel.centerOfGeometry()

        for residue in self.residue_atoms :
            com2 = xyz[residue.first:residue.last+1].mean(axis=0)
            dr = pbc.unwrap_vector(com1 - com2, self.processor.currbox)
            xyz[residue.first:residue.last+1] = xyz[residue.first:residue.last+1] + dr

        delta = com1 - self.processor.currbox/2.0
        self.processor.currsnap._pos = xyz - delta

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
        parser.add_argument('--usedict',action="store_true",help="if to use dictionary vectors",default=False)

    def setup(self,args):
        protsel = self.processor.universe.selectAtoms(args.pmask)
        self.usedict = args.usedict
        if not args.usedict:
            self.atm2 = protsel.selectAtoms("name "+args.atoms[1])
            self.atm1 = protsel.selectAtoms("name "+args.atoms[0]+
                                        " and byres name "+args.atoms[1])
        else:
            atm1 = []
            atm2 = []
            for res in self.processor.universe.selectAtoms("protein").residues:
                if res.name == 'VAL' :
                  atm1.append(res.CB)
                  atm1.append(res.CB)
                  atm2.append(res.CG1)
                  atm2.append(res.CG2)
                elif res.name == 'SER' :
                  atm1.append(res.CB)
                  atm1.append(res.CB)
                  atm2.append(res.HB2)
                  atm2.append(res.HB3)
                elif res.name == 'THR' :
                  atm1.append(res.CB)
                  atm2.append(res.CG2)
                elif res.name == 'ILE' :
                  atm1.append(res.CG1)
                  atm2.append(res.CD1)
                elif res.name == 'LEU' :
                  atm1.append(res.CG)
                  atm1.append(res.CG)
                  atm2.append(res.CD1)
                  atm2.append(res.CD2)
                elif res.name == 'MET' :
                  atm1.append(res.SD)
                  atm2.append(res.CE)
                elif res.name == 'ASN' :
                  atm1.append(res.ND2)
                  atm1.append(res.ND2)
                  atm2.append(res.HD21)
                  atm2.append(res.HD22)
                elif res.name == 'GLN' :
                  atm1.append(res.NE2)
                  atm1.append(res.NE2)
                  atm2.append(res.HE21)
                  atm2.append(res.HE22)
                elif res.name == 'PHE' :
                  atm1.append(res.CD1)
                  atm2.append(res.HD1)
                elif res.name == 'HID' or res.name == 'HIE' or res.name == 'HIP' :
                  atm1.append(res.CD2)
                  atm2.append(res.HD2)
                elif res.name == 'TYR' :
                  atm1.append(res.CD1)
                  atm2.append(res.HD1)
                elif res.name == 'PRO' :
                  atm1.append(res.CG)
                  atm1.append(res.CG)
                  atm2.append(res.HG2)
                  atm2.append(res.HG3)
                elif res.name == 'LYS' :
                  atm1.append(res.CB)
                  atm1.append(res.CB)
                  atm1.append(res.CG)
                  atm1.append(res.CG)
                  atm1.append(res.CD)
                  atm1.append(res.CD)
                  atm1.append(res.CE)
                  atm1.append(res.CE)
                  atm2.append(res.HB2)
                  atm2.append(res.HB3)
                  atm2.append(res.HG2)
                  atm2.append(res.HG3)
                  atm2.append(res.HD2)
                  atm2.append(res.HD3)
                  atm2.append(res.HE2)
                  atm2.append(res.HE3)
                elif res.name == 'ARG' :
                  atm1.append(res.CB)
                  atm1.append(res.CB)
                  atm1.append(res.CG)
                  atm1.append(res.CG)
                  atm1.append(res.CD)
                  atm1.append(res.CD)
                  atm1.append(res.NE)
                  atm2.append(res.HB2)
                  atm2.append(res.HB3)
                  atm2.append(res.HG2)
                  atm2.append(res.HG3)
                  atm2.append(res.HD2)
                  atm2.append(res.HD3)
                  atm2.append(res.HE)
                elif res.name == 'ASP' :
                  atm1.append(res.CB)
                  atm1.append(res.CB)
                  atm2.append(res.HB2)
                  atm2.append(res.HB3)
                elif res.name == 'GLU' :
                  atm1.append(res.CG)
                  atm1.append(res.CG)
                  atm2.append(res.HG2)
                  atm2.append(res.HG3)
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
                if not self.usedict:
                    f.write("%s %d %.5f\n"%(atm.resname,atm.resnum+self.resoffset,rs2.mean()))
                else :
                    if reslist and atm.resnum != prevres:
                        av = np.asarray(reslist).mean()
                        f.write("%s %d %.5f\n"%(prevatom.resname,prevatom.resnum+self.resoffset,av))
                        reslist = []
                    reslist.append(rs2.mean())
                    prevres = atm.resnum
                    prevatom = atm
            if self.usedict:
                av = np.asarray(reslist).mean()
                f.write("%s %d %.5f\n"%(atm.resname,atm.resnum+self.resoffset,av))

class MempropAnalysis(TrajectoryAction):
    def add_arguments(self, parser):
        parser.add_argument('--pmask',help="the selectiom mask for phosphor atoms",default="name P")
        parser.add_argument('--lipidmask',help="the selectiom mask for lipid residues",default="resname POPC")
        parser.add_argument('--watmask',help="the selectiom mask for water residues",default="resname SOL")
        parser.add_argument('--watvol',type=float,help="the volume of a water molecule in nm3",default=0.0306)
        parser.add_argument('-o','--out',help="the output prefix",default="memprop")

    def setup(self,args):
        self.out = args.out
        self.phosphorsel = self.processor.universe.selectAtoms(args.pmask)
        self.lipidsel = self.processor.universe.selectAtoms(args.lipidmask)
        watsel  = self.processor.universe.selectAtoms(args.watmask)
        self.nlipid = len(self.lipidsel.residues)
        self.nwat = len(watsel.residues)
        self.watvol = args.watvol
        nphosph = len(self.phosphorsel.residues)
        if self.nlipid == 0 or self.nwat == 0 or nphosph == 0 :
            raise Exception("Either number of lipids (%d), water (%d) or phosphor atoms (%d)  is zero"%(self.nlipid,self.nwat,nphosph))
        self.apllist = []
        self.vpllist = []
        # Setup edges to cover the entire simulation box
        zpos = self.processor.universe.coord._pos[:,2] - self.lipidsel.positions[:,2].mean()
        self.resolution = 0.25
        self.edges = np.arange(zpos.min(),zpos.max()+self.resolution,self.resolution)
        self.density = np.zeros(self.edges.shape)
        self.records = []

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
        self.records.append(MDRecord(self.processor.currtime,[self.apllist[-1],self.vpllist[-1],self._calc_dhh()]))

    def finalize(self):
        """
        Calculate average APL and VPL as well as distance
        between peaks in the phosphor density
        """
        dhh = self._calc_dhh()
        apl = np.asarray(self.apllist).mean()
        vpl = np.asarray(self.vpllist).mean()
        with open(self.out+".txt","w") as f :
            f.write("%.3f\t%.3f\t%.3f\n"%(apl,vpl,dhh))
        self._write_records(postfix="_dt.txt")

    def _calc_dhh(self) :
        mid = int(self.density.shape[0]/2)
        dens_first = self.density[:mid]
        dens_last = self.density[mid:]
        max_first = np.argmax(dens_first)
        max_last = np.argmax(dens_last)
        return (max_last + mid - max_first) / 10.0 * self.resolution

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
        self.selection = self.processor.universe.selectAtoms(args.mask)
        self.masses = np.asarray([atom.mass for atom in self.selection])
        self.normal = np.asarray(args.normal)
        self.records = []
        self.out = args.out

    def process(self):
        xyz = pbc.make_whole_xyz(self.selection.get_positions(),
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
        parser.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
        parser.add_argument('-o','--out',help="the output",default="rmsf.txt")

    def setup(self, args):
        self.refuni = md.Universe(self.processor.args.struct)
        self.protsel = self.processor.universe.selectAtoms(args.pmask)
        self.atoms = args.atoms
        self.alignmask = "%s"%(" or ".join("name %s"%a for a in args.atoms))
        natm = len(self.processor.universe.atoms)
        self.sumcoords2 = np.zeros([natm,3])
        self.sumcoords = np.zeros([natm,3])
        self.out = args.out

    def process(self):
        rmsd = align.alignto(self.processor.universe, self.refuni, select=self.alignmask)
        xyz = self.processor.currsnap._pos
        self.sumcoords += xyz
        self.sumcoords2 += xyz*xyz

    def finalize(self):
        nsnap = float(self.processor.currtime / self.processor.dt)
        self.sumcoords = self.sumcoords / nsnap
        self.sumcoords2 = self.sumcoords2 / nsnap
        var = self.sumcoords2 - (self.sumcoords * self.sumcoords)
        bfac = (8.0/3.0)*np.pi*np.pi*var.sum(axis=1)

        with open(self.out,'w') as f:
            for residue in self.protsel.residues:
                rbfac = 0.0
                masssum = 0.0
                for atom in residue:
                    if atom.name not in self.atoms : continue
                    rbfac = bfac[atom.number]*atom.mass
                    masssum += atom.mass
                f.write("%d %.3f\n"%(residue.id,rbfac / masssum))

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
        self.atom1 = self.processor.universe.selectAtoms(args.atoms[0])
        if len(self.atom1) == 0:
            raise Exception("Selection %s did not select any atoms"%args.atoms[0])
        self.atom2 = self.processor.universe.selectAtoms(args.atoms[1])
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
