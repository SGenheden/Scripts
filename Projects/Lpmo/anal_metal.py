# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to analyse the metal site in LPMO simulations

Examples:

"""

import numpy as np
import MDAnalysis as md
import MDAnalysis.lib.distances as mddist
import MDAnalysis.analysis.align as align

from sgenlib import moldyn
from sgenlib import mdactions
from sgenlib import geo

class MetalSiteAnalysis(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--copper', help="the selection mask for the copper", default="name CU")
        parser.add_argument('--bonds', nargs="+", help="the selection mask for the bonds to consider",
            default=["name N and resid 0", "name ND1 and resid 0",
                     "name NE2 and resid 77", "name OW and resid 236",
                     "name OH and resid 163"])
        parser.add_argument('--watsel', nargs="+", default=["name OW and resid 236", "name OW and resid 237"])
        parser.add_argument('--coordsel', nargs=4, help="the selection mask for the water coordination",
                            default=["name N and resid 0", "name ND1 and resid 0", "name NE2 and resid 77"])
        parser.add_argument('--water', help="the selection mask for the water", default="name OW and (resid 237-8072)")
        parser.add_argument('-o','--out',help="the output",default="metal_anal")

    def setup(self, args):
        self.copper = self.processor.universe.select_atoms(args.copper)
        self.bonds = []
        self.bonddata = []
        for bond in args.bonds:
            sel = self.processor.universe.select_atoms("%s or (%s)"%(args.copper,bond))
            print bond, len(sel)
            if len(sel) == 2 :
                self.bonds.append(sel)
        self.water = self.processor.universe.select_atoms(args.water)
        print "Selected %d water atoms"%len(self.water)
        self.incontact = []
        self.coordsel = self.processor.universe.select_atoms(" or ".join(["(%s)"%s for s in args.coordsel]))
        print "Coordination selection contains %d atoms"%len(self.coordsel)
        self.coordsel_wat = [self.processor.universe.select_atoms(s) for s in args.watsel]
        self.coordination = []
        self.wcoordination = []
        self.cutoff = 4.0
        self.rmsds = []

        self.out = args.out

    def process(self):

        def _calc_dist(bond):
            r = bond[0].pos - bond[1].pos
            return np.sqrt(np.sum(r * r))*0.1

        def  _calc_coord(xyz, normal, origin):
            return geo.angle(normal, origin - xyz)*180.0/np.pi
        self.bonddata.append(mdactions.MDRecord(self.processor.currtime,
                                [_calc_dist(b) for b in self.bonds]))

        pos = self.coordsel.positions
        normal = np.cross(pos[1,:]-pos[0,:],pos[2,:]-pos[0,:])
        origin = pos.mean(axis=0)
        self.coordination.append(mdactions.MDRecord(self.processor.currtime,
                    [_calc_coord(w.positions, normal, origin) for w in self.coordsel_wat]))

        r = self.water.positions - self.copper.positions
        d = sel = np.sqrt(np.sum(r * r, axis=1))
        self.wcoordination.extend([_calc_coord(self.water.positions[i,:], normal, origin) \
                    for i, d in enumerate(d) if d <= self.cutoff])
        n = np.sum(d <= self.cutoff)
        self.incontact.append(mdactions.MDRecord(self.processor.currtime,n))


    def finalize(self):
        self.records = self.bonddata
        self._write_records(postfix="_bonds.txt")

        self.records = self.incontact
        self._write_records(postfix="_wcontacts.txt")

        self.records = self.coordination
        self._write_records(postfix="_coordination.txt")

        self.wcoordination = np.asarray(self.wcoordination)
        with open(self.out+"_wcoordination.txt", "w") as f :
            for i, c in enumerate(self.wcoordination):
                f.write("%d %.3f\n"%(i+1, c))

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Analyse metal site")
    analysis = MetalSiteAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
