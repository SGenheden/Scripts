# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate APL, VPL and Dhh from an MD trajectory

Prints out the results in SI units

Examples:
    md_memprop.py -f md2.xtc -s md2.gro
    md_memprop.py -f md2.xtc -s md2.gro --watmask "resname W" --pmask "name PO4" --watvol 0.12
"""

import sys

import numpy as np

from sgenlib import moldyn

class MempropAnalysis(moldyn.AnalysisAction):
    def __init__(self,processor):
        super(MempropAnalysis,self).__init__(processor)
        self.outname = processor.args.out
        self.phosphorsel = processor.universe.selectAtoms(processor.args.pmask)
        self.lipidsel = processor.universe.selectAtoms(processor.args.lipidmask)
        watsel  = processor.universe.selectAtoms(processor.args.watmask)
        self.nlipid = len(self.lipidsel.residues)
        self.nwat = len(watsel.residues)
        self.watvol = processor.args.watvol
        nphosph = len(self.phosphorsel.residues)
        if self.nlipid == 0 or self.nwat == 0 or nphosph == 0 :
            raise Exception("Either number of lipids (%d), water (%d) or phosphor atoms (%d)  is zero"%(self.nlipid,self.nwat,nphosph))
        self.apllist = []
        self.vpllist = []
        # Setup edges to cover the entire simulation box
        zpos = processor.universe.coord._pos[:,2] - self.lipidsel.positions[:,2].mean()
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

if __name__ == '__main__' :

    #print " ".join(sys.argv)
    processor = moldyn.TrajectoryProcessor("Calculate membrane properties")
    processor.argparser.add_argument('--pmask',help="the selectiom mask for phosphor atoms",default="name P")
    processor.argparser.add_argument('--lipidmask',help="the selectiom mask for lipid residues",default="resname POPC")
    processor.argparser.add_argument('--watmask',help="the selectiom mask for water residues",default="resname SOL")
    processor.argparser.add_argument('--watvol',type=float,help="the volume of a water molecule in nm3",default=0.0306)
    processor.argparser.add_argument('-o','--out',help="the output filename",default="memprop.txt")
    processor.setup()

    analysis = MempropAnalysis(processor)
    processor.process()
