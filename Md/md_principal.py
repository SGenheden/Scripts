# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate angle between the principal axis
of a selection and an arbitrary vector

By default is calculates the principal axis of the CA atoms
towards the z-axis.

Examples:
    md_principal.py -f sim.dcd -s ref.pdb
"""

import sys

import numpy as np

from sgenlib import pbc
from sgenlib import geo
from sgenlib import moldyn

class PrincipalAxisAnalysis(moldyn.AnalysisAction):
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
    def __init__(self,processor):
        super(PrincipalAxisAnalysis,self).__init__(processor)
        self.selection = processor.universe.selectAtoms(processor.args.mask)
        self.masses = np.asarray([atom.mass for atom in self.selection])
        self.normal = np.asarray(processor.args.normal)
        self.records = []

    def process(self):
        xyz = pbc.make_whole_xyz(self.selection.get_positions(),
                                    self.processor.currbox)
        moi = geo.moment_of_inertia(xyz-xyz.mean(axis=0),self.masses)
        princip = geo.principal_axes(moi)
        alpha = geo.angle(princip[0,:],self.normal)
        dalpha = pbc.unwrap_vector(alpha,np.pi)
        alpha = np.abs(alpha-dalpha)*180.0/np.pi
        self.records.append(moldyn.MDRecord(self.processor.currtime,alpha))

    def finalize(self):
        """
        Write out average alpha and then all alphas to disc
        """
        alphas = np.asarray([entry.value for entry in self.records])
        print "Mean = %.3f Std = %.3f"%(alphas.mean(),alphas.std())
        self._write_records()

if __name__ == '__main__' :

    print " ".join(sys.argv)
    processor = moldyn.TrajectoryProcessor("Calculte the principcal axis of a molecule and it angles with an axis")
    processor.argparser.add_argument('-m','--mask',help="the selectiom mask",default="name CA")
    processor.argparser.add_argument('-n','--normal',type=float,nargs=3,help="the normal vector",default=[0.0,0.0,1.0])
    processor.argparser.add_argument('-o','--out',help="the output filename",default="alpha.txt")
    processor.setup()

    analysis = PrincipalAxisAnalysis(processor)
    processor.process()
