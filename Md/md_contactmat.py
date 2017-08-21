# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to compute the contact matrix

Examples:

"""

import numpy as np
import MDAnalysis as md
import MDAnalysis.analysis.distances as mddist
from MDAnalysis.core.groups import AtomGroup

from sgenlib import moldyn
from sgenlib import mdactions

class ContactMatrixAnalysis(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('-o','--out',help="the output",default="contactmat.txt")

    def setup(self, args):
        self.calphas = self.processor.universe.select_atoms("name CA")
        indices = []
        for res in self.calphas.residues :
            try :
                indices.append(res["CB"].index)
            except :
                indices.append(res["CA"].index)
        self.cbetas = AtomGroup(indices, self.processor.universe)
        n = len(self.calphas)
        self.contacts = np.zeros(n*(n-1)/2)
        self.tempmat = np.zeros(n*(n-1)/2)
        self.out = args.out


    def process(self):

        proj = self.calphas.positions +  \
            1.5*(self.cbetas.positions-self.calphas.positions)
        mddist.self_distance_array(proj, box=self.processor.currbox,
                                        result=self.tempmat)
        self.contacts += self.tempmat

    def finalize(self):
        self.contacts /= float(self.processor.nprocessed)
        self.contacts *= 0.1 # To  Nm
        outmat = np.zeros([len(self.calphas),len(self.calphas)])
        k = 0
        for i in range(len(self.calphas)) :
            for j in range(i+1, len(self.calphas)) :
                outmat[i,j] = self.contacts[k]
                outmat[j,i] = self.contacts[k]
                k += 1

        with open(self.out, "w") as f :
            for i in range(outmat.shape[0]) :
                for j in range(outmat.shape[1]) :
                    f.write("%.3f\t"%outmat[i,j])
                f.write("\n")

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate contact matrix")
    analysis = ContactMatrixAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
