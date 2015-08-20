# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to transform a molecular dynamics trajectory.

It can
1) Make molecules whole over periodic boxes
2) Center a protein in the middle of the box
3) Align protein with an RMSD fit

Works only for rectangular boxes.

Examples:
    md_centerwhole.py -f sim.dcd -s ref.pdb
    md_centerwhole.py -f sim.dcd -s ref.pdb --noalign
"""

import sys
from collections import namedtuple

import MDAnalysis as md
import MDAnalysis.analysis.align as align
import numpy as np

from sgenlib import pbc
from sgenlib import moldyn

ResidueAtoms = namedtuple("ResidueAtoms",["first","last"])

class CenterWholeAlign(moldyn.TrajectoryAction):
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
    records : list of moldyn.MDRecords
        the RMSD record at each processes snapshot
    writer : MDAnalysis.Writer
        the output trajectory writer
    """
    def __init__(self,processor):
        super(CenterWholeAlign,self).__init__(processor)
        self.refuni = md.Universe(processor.args.struct)
        self.protsel = processor.universe.selectAtoms(processor.args.pmask)
        if len(self.protsel) == 0 :
            self.processor.args.nocenter = True
            self.processor.args.noalign = True

        self.residues = []
        self.residue_atoms = []
        for res in processor.universe.selectAtoms("not "+processor.args.pmask).residues:
            if len(res) > 1 :
                self.residues.append(res)
            self.residue_atoms.append(ResidueAtoms(res[0].number,res[-1].number))

        self.records = []
        self.writer = md.Writer(processor.args.out+".dcd",
                                    processor.universe.trajectory.numatoms)

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
            self.records.append(moldyn.MDRecord(self.processor.currtime,rmsd))

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

    @classmethod
    def add_arguments(cls,processor):
        processor.argparser.add_argument('--bbmask',help="the selectiom mask for backbone",default="name CA")
        processor.argparser.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
        processor.argparser.add_argument('-o','--out',help="the output",default="centerwhole")
        processor.argparser.add_argument('--noalign',action="store_true",help="turns off alignment",default=False)
        processor.argparser.add_argument('--nocenter',action="store_true",help="turns off centering",default=False)

if __name__ == '__main__' :

  processor = moldyn.TrajectoryProcessor("Center and make a trajectory whole")
  CenterWholeAlign.add_arguments(processor)
  processor.setup(printargs=True)

  analysis = CenterWholeAlign(processor)
  processor.process()
