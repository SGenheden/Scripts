
from collections import namedtuple

import MDAnalysis as md
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
    def __init__(self,processor):
        self.processor = processor
        self.dosubsample = False
        processor.actions.append(self)

    def add_arguments(self):
        """
        Function that is called to add action-specific arguments
        """
        pass

    def setup(self):
        """
        Function that is called after the processor has parsed the command-line
        arguments.
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
            with open(self.processor.args.out+postfix,'w') as f :
                for entry in self.records:
                    f.write("%.0f %.3f\n"%(entry.time,entry.value))


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
    def add_arguments(self):
        self.processor.add_argument('--bbmask',help="the selectiom mask for backbone",default="name CA")
        self.processor.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
        self.processor.add_argument('-o','--out',help="the output",default="centerwhole")
        self.processor.add_argument('--noalign',action="store_true",help="turns off alignment",default=False)
        self.processor.add_argument('--nocenter',action="store_true",help="turns off centering",default=False)

    def setup(self):
        self.refuni = md.Universe(self.processor.args.struct)
        self.protsel = self.processor.universe.selectAtoms(self.processor.args.pmask)
        if len(self.protsel) == 0 :
            self.processor.args.nocenter = True
            self.processor.args.noalign = True

        self.residues = []
        self.residue_atoms = []
        for res in self.processor.universe.selectAtoms("not "+self.processor.args.pmask).residues:
            if len(res) > 1 :
                self.residues.append(res)
            self.residue_atoms.append(ResidueAtoms(res[0].number,res[-1].number))

        self.records = []
        self.writer = md.Writer(self.processor.args.out+".dcd",
                                    self.processor.universe.trajectory.numatoms)

    def process(self):
        if len(self.protsel) > 0:
            xyz = pbc.make_whole_xyz(self.protsel.get_positions(),self.processor.currbox)
            self.protsel.set_positions(xyz)
        for res in self.residues :
            xyz = pbc.make_whole_xyz(res.get_positions(),self.processor.currbox)
            res.set_positions(xyz)
        if not self.processor.args.nocenter :
            self._center()
        if not self.processor.args.noalign :
            rmsd = align.alignto(self.processor.universe, self.refuni,
                                select=self.processor.args.bbmask)[1]
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
        com1 = xyz[self.protsel[0].number:self.protsel[-1].number+1].mean(axis=0)

        for residue in self.residue_atoms :
            com2 = xyz[residue.first:residue.last+1].mean(axis=0)
            dr = pbc.unwrap_vector(com1 - com2, self.processor.currbox)
            xyz[residue.first:residue.last+1] = xyz[residue.first:residue.last+1] + dr

        delta = com1 - self.processor.currbox/2.0
        xyz = xyz - delta

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
    def add_arguments(self):
        self.processor.add_argument('--atoms',nargs=2,help="the atom names making the vectors",default=["N","H"])
        self.processor.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
        self.processor.add_argument('-o','--out',help="the output name",default="s2.txt")

    def setup(self):
        protsel = self.processor.universe.selectAtoms(self.processor.args.pmask)
        self.atm2 = protsel.selectAtoms("name "+self.processor.args.atoms[1])
        self.atm1 = protsel.selectAtoms("name "+self.processor.args.atoms[0]+
                                        " and byres name "+self.processor.args.atoms[1])
        self.mat = np.zeros([len(self.atm1),len(self.atm1)])
        self.s2list = []
        self.outname = self.processor.args.out
        self.dosubsample = True

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
            for i,(atm,rs2) in enumerate(zip(self.atm1,self.s2list.T)):
                f.write("%d %.5f\n"%(i,rs2.mean()))

class MempropAnalysis(TrajectoryAction):
    def add_arguments(self):
        self.processor.add_argument('--pmask',help="the selectiom mask for phosphor atoms",default="name P")
        self.processor.add_argument('--lipidmask',help="the selectiom mask for lipid residues",default="resname POPC")
        self.processor.add_argument('--watmask',help="the selectiom mask for water residues",default="resname SOL")
        self.processor.add_argument('--watvol',type=float,help="the volume of a water molecule in nm3",default=0.0306)
        self.processor.add_argument('-o','--out',help="the output filename",default="memprop.txt")

    def setup(self):
        self.outname = self.processor.args.out
        self.phosphorsel = self.processor.universe.selectAtoms(self.processor.args.pmask)
        self.lipidsel = self.processor.universe.selectAtoms(self.processor.args.lipidmask)
        watsel  = self.processor.universe.selectAtoms(self.processor.args.watmask)
        self.nlipid = len(self.lipidsel.residues)
        self.nwat = len(watsel.residues)
        self.watvol = self.processor.args.watvol
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

    def finalize(self):
        """
        Calculate average APL and VPL as well as distance
        between peaks in the phosphor density
        """
        mid = int(self.density.shape[0]/2)
        dens_first = self.density[:mid]
        dens_last = self.density[mid:]
        max_first = np.argmax(dens_first)
        max_last = np.argmax(dens_last)
        dhh = (max_last + mid - max_first) / 10.0 * self.resolution
        apl = np.asarray(self.apllist).mean()
        vpl = np.asarray(self.vpllist).mean()
        with open(self.outname,"w") as f :
            f.write("%.3f\t%.3f\t%.3f\n"%(apl,vpl,dhh))


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
    def add_arguments(self):
        self.processor.add_argument('-m','--mask',help="the selectiom mask",default="name CA")
        self.processor.add_argument('-n','--normal',type=float,nargs=3,help="the normal vector",default=[0.0,0.0,1.0])
        self.processor.add_argument('-o','--out',help="the output filename",default="alpha.txt")

    def setup(self):
        self.selection = self.processor.universe.selectAtoms(self.processor.args.mask)
        self.masses = np.asarray([atom.mass for atom in self.selection])
        self.normal = np.asarray(self.processor.args.normal)
        self.records = []

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
    def add_arguments(self):
        self.processor.add_argument('--atoms',nargs="+",help="the atom names in the  backbone",default=["CA","N","C"])
        self.processor.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
        self.processor.add_argument('-o','--out',help="the output",default="rmsf.txt")

    def setup(self):
        self.refuni = md.Universe(self.processor.args.struct)
        self.protsel = self.processor.universe.selectAtoms(self.processor.args.pmask)
        self.alignmask = "%s"%(" or ".join("name %s"%a for a in self.processor.args.atoms))
        natm = len(self.processor.universe.atoms)
        self.sumcoords2 = np.zeros([natm,3])
        self.sumcoords = np.zeros([natm,3])

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

        with open(self.processor.args.out,'w') as f:
            for residue in self.protsel.residues:
                rbfac = 0.0
                masssum = 0.0
                for atom in residue:
                    if atom.name not in self.processor.args.atoms : continue
                    rbfac = bfac[atom.number]*atom.mass
                    masssum += atom.mass
                f.write("%d %.3f\n"%(residue.id,rbfac / masssum))
