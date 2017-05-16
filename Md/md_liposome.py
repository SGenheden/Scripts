# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the radius and area per lipid of a liposome

Examples:
  md_liposome.py -s ref.gro -f sim.xtc --inner "name PO4 and resid 8161:9107" --outer "name PO4 and resid 9108:10688"
"""

import numpy as np
import MDAnalysis as md
import MDAnalysis.lib.distances as mddist

from sgenlib import moldyn
from sgenlib import mdactions
from sgenlib import pbc
from sgenlib import geo

class LiposomeAnalysis(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--inner', help="the selectiom mask for the inner lipids")
        parser.add_argument('--outer', help="the selectiom mask for the outer lipids")
        parser.add_argument('--lipids', help="the selection mask for all lipids")
        parser.add_argument('-o','--out',help="the output",default="lipoanal.txt")
        parser.add_argument('-of','--outfile',help="the output trajectory")

    def setup(self, args):
        self.inner = self.processor.universe.select_atoms(args.inner)
        print "Inner selection contains %d particles"%len(self.inner.residues)
        self.outer = self.processor.universe.select_atoms(args.outer)
        print "Outer selection contains %d particles"%len(self.outer.residues)
        self.out = args.out
        self.records = []
        self.innerfac = 4*np.pi/float(len(self.inner.residues))
        self.outerfac = 4*np.pi/float(len(self.outer.residues))
        self.lipids = self.processor.universe.select_atoms(args.lipids)
        print "Total number of lipids: %d"%len(self.lipids.residues)
        if args.outfile is not None :
            self.writer = md.Writer(args.outfile,
                                    self.processor.universe.trajectory.n_atoms)
            self.residue_atoms = []
            for res in self.processor.universe.select_atoms("all").residues:
                self.residue_atoms.append(mdactions.ResidueAtoms(res[0].number,res[-1].number))
        else:
            self.writer = None


    def process(self):

        if self.writer is not None:
            # Make liposome whole
            residues = self.lipids.residues
            for res2 in residues :
                dr = residues[0].center_of_geometry()-res2.center_of_geometry()
                dr = pbc.unwrap_vector(dr, self.processor.currbox)
                res2.set_positions(res2.get_positions()+dr)

            # Center liposome
            xyz = self.processor.currsnap._pos
            com1 = np.median(self.lipids.get_positions(),axis=0)
            for residue in self.residue_atoms :
                com2 = xyz[residue.first:residue.last+1].mean(axis=0)
                dr = pbc.unwrap_vector(com1 - com2, self.processor.currbox)
                xyz[residue.first:residue.last+1] = xyz[residue.first:residue.last+1] + dr
            delta = com1 - self.processor.currbox/2.0
            self.processor.currsnap._pos = xyz - delta

            self.writer.write(self.processor.currsnap)

        innercom = self.inner.center_of_geometry()
        outercom = self.outer.center_of_geometry()
        lipcom = np.asarray([0.5*(innercom+outercom)])
        rinner = mddist.distance_array(lipcom,self.inner.positions,None).mean()*0.1
        router = mddist.distance_array(lipcom,self.outer.positions,None).mean()*0.1
        eigval, asphericity = geo.sphericity(self.inner.positions)
        axdevinner = np.sqrt(eigval.max()/eigval.mean())*100
        eigval, asphericity = geo.sphericity(self.outer.positions)
        axdevout = np.sqrt(eigval.max()/eigval.mean())*100
        self.records.append(mdactions.MDRecord(self.processor.currtime,
            [rinner,router,rinner**2*self.innerfac,router**2*self.outerfac, axdevinner, axdevout]))

    def finalize(self):
        self._write_records()

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate liposome properties")
    analysis = LiposomeAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
