# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate the RDF of a group of atoms

Experimental code!

"""

import numpy as np
import MDAnalysis.core.distances as mddist

from sgenlib import moldyn
from sgenlib import mdactions
from sgenlib import pbc

class RdfAnalysis(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--ref',help="the selectiom mask for the central group")
        parser.add_argument('--sels',nargs="+",help="the selectiom mask for the selection group")
        parser.add_argument('--rmax',type=float,help="the maximum distance")
        parser.add_argument('--usecom',action="store_true",help="calculate RDF between com",default=False)
        parser.add_argument('-o','--out',help="the output",default="rdf.txt")

    def setup(self, args):
        self.ref = self.processor.universe.selectAtoms(args.ref)
        print "Reference selection contains %d atoms"%len(self.ref)
        self.sels = [self.processor.universe.selectAtoms(s) for s in args.sels]
        self.out = args.out
        if args.rmax is not None:
            self.rmax = args.rmax
        else:
            self.rmax = self.processor.universe.coord.dimensions[:2].max()
        self.resolution = 0.1
        self.edges = np.arange(0.0,self.rmax+self.resolution,self.resolution)
        self.densities = [np.zeros(self.edges.shape[0]-1) for s in self.sels]
        if not args.usecom:
            self.dists = [np.zeros((len(self.ref),len(s))) for s in self.sels]
        self.usecom = args.usecom
        self.vol = 0.0

    def process(self):

        self.vol += self.processor.currsnap.volume
        if self.usecom:
            refcom = np.asarray([self.ref.centerOfGeometry()])
            for i,sel in enumerate(self.sels):
                selcom = np.asarray([r.centerOfGeometry() for r in sel.residues])
                dist = mddist.distance_array(refcom, sel.positions, self.processor.currbox)
                h, e = np.histogram(dist,bins=self.edges)
                self.densities[i] += h
        else:
            for i,(sel,dist) in enumerate(zip(self.sels,self.dists)):
                mddist.distance_array(self.ref.positions,sel.positions,
                        self.processor.currbox,dist)
                h, e = np.histogram(dist,bins=self.edges)
                self.densities[i] += h

    def finalize(self):

        boxvol = 1.0 / self.vol
        print self.vol / self.processor.nprocessed
        radii = 0.5 * (self.edges[1:] + self.edges[:-1])
        vol = 4.0 / 3.0 * np.pi * (np.power(self.edges[1:], 3) - np.power(self.edges[:-1], 3))
        if self.usecom:
            norms = [(len(s.residues)) * boxvol * vol for s in self.sels]
        else:
            norms = [(len(s)*len(self.ref)) * boxvol * vol for s in self.sels]
        print norms
        rdfs = [den / norm  for den,norm in zip(self.densities,norms)]
        rdfs = np.asarray(rdfs)
        with open(self.out,'w') as f :
            for i,rdfr in enumerate(rdfs.T):
                f.write("%.2f\t%s\n"%(i*self.resolution,"\t".join(["%.3f"%r for r in rdfr])))
        #self._write_records(postfix=".hist")

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate the RDF of a group of atoms")
    analysis = RdfAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
