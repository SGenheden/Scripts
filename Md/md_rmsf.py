# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate RMSF of an MD trajectory

By defaults calculates the RMSF of the C, CA and N atoms

Examples:
    md_rmsf.py -f sim.dcd -s ref.pdb
"""

import sys

import MDAnalysis as md
import MDAnalysis.analysis.align as align
import numpy as np

from sgenlib import moldyn

class RMSFAnalysis(moldyn.TrajectoryAction):
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
    def __init__(self,processor):
        super(RMSFAnalysis,self).__init__(processor)
        self.refuni = md.Universe(processor.args.struct)
        self.protsel = processor.universe.selectAtoms(processor.args.pmask)
        self.alignmask = "%s"%(" or ".join("name %s"%a for a in processor.args.atoms))
        natm = len(processor.universe.atoms)
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

    @classmethod
    def add_arguments(cls,processor):
        processor.argparser.add_argument('--atoms',nargs="+",help="the atom names in the  backbone",default=["CA","N","C"])
        processor.argparser.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
        processor.argparser.add_argument('-o','--out',help="the output",default="rmsf.txt")

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate the bfactor of atoms")
    RMSFAnalysis.add_arguments(processor)
    processor.setup(printargs=True)

    analysis = RMSFAnalysis(processor)
    processor.process()
