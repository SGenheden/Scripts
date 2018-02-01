# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to perform MSD analysis of a trajectory

Experimental code!
"""

import numpy as np

from sgenlib import moldyn
from sgenlib import mdactions
from sgenlib import pbc

class MsdWalker(object):

    def __init__(self, time, owner):
        self.time = time
        self.origin = owner.atoms.get_positions()
        self.rtrue = np.array(self.origin,copy=True)
        self.processor = owner.processor
        self.owner = owner
        self.inplane = owner.inplane
        self.direction = owner.direction

    def update(self):
        if (self.processor.currtime - self.time) % self.processor.freq == 0 :
            positions = self.owner.atoms.get_positions()
            idx = int((self.processor.currtime - self.time) / self.processor.freq - 1)
            #print self.processor.currtime, self.time, idx
            dr = pbc.unwrap_vector(self.rtrue - positions, self.processor.currbox)
            self.rtrue = positions+dr
            dr = self.rtrue - self.origin
            if self.owner.varmask:
                dr = dr[self.owner.mask,:]
            if not self.inplane and not self.direction :
                self.owner.msds[idx].value.append(np.mean((dr**2).sum(axis=1)))
            elif self.inplane:
                self.owner.msds[idx].value.append(np.mean((dr[:,:2]**2).sum(axis=1)))
            elif self.direction:
                self.owner.msds[idx].value.append(np.mean((dr[:,2]**2)))

class MsdAnalysis(mdactions.TrajectoryAction) :

    def add_arguments(self, parser):
        parser.add_argument('--mask',help="the selectiom mask for the group")
        parser.add_argument('--wfreq',type=int,help="the walker frequency")
        parser.add_argument('-o','--out',help="the output",default="msd.txt")
        parser.add_argument('--diff',action="store_true",help="turn on estimate of D",default=False)
        parser.add_argument('--variable',action="store_true",help="turn on variable mask",default=False)
        parser.add_argument('--inplane',action="store_true",help="turn on in plane calculation",default=False)
        parser.add_argument('--direction',action="store_true",help="turn on in direction calculation",default=False)

    def setup(self, args):
        self.varmask = args.variable
        if self.varmask :
            self.maskstr = args.mask
            self.atoms = self.processor.universe.select_atoms("all")
        else:
            self.atoms = self.processor.universe.select_atoms(args.mask)
            print "%d atoms selected"%len(self.atoms)
        self.out = args.out
        self.records = []
        self.dosubsample = True
        self.wfreq = -1 if args.wfreq is None else args.wfreq
        self.walkers = []
        self.msds = [] #[[]]
        self.dodiff = args.diff
        self.inplane = args.inplane
        self.direction = args.direction

    def process(self):

        if self.processor.currtime % self.processor.freq == 0:
            self.msds.append(mdactions.MDRecord(self.processor.currtime,[]))

        if self.varmask:
            self.selatoms = self.processor.universe.select_atoms(self.maskstr)
            self.mask = np.zeros(len(self.atoms),dtype=bool)
            for atom in self.selatoms:
                self.mask[atom.number] = True
            self.records.append(mdactions.MDRecord(self.processor.currtime,np.sum(self.mask)))

        if len(self.msds) > 0 :
            for walker in self.walkers:
                walker.update()

        if self.processor.currtime == self.processor.dt or \
            (self.wfreq > 0 and self.processor.currtime % self.wfreq == 0):
            self.walkers.append(MsdWalker(self.processor.currtime -
                                            self.processor.dt,self))

    def subsample(self):
        #if len(self.msds) == 1:
        #    print "here: ", self.processor.currtime, self.msds[0]
        #self.msds.append(mdactions.MDRecord(self.processor.currtime,[]))
        pass

    def finalize(self):
        self._write_records(".natm")
        self.records = []
        for msd  in self.msds[1:] :
            av = np.asarray(msd.value).mean()
            std = np.asarray(msd.value).std() / np.sqrt(len(msd))
            self.records.append(mdactions.MDRecord(msd.time, [av,std]))
        if self.dodiff :
            msd = np.asarray([e.value[0] for e in self.records])
            time = np.asarray([e.time for e in self.records])
            imin = int(np.floor(0.1*msd.shape[0]))
            imax = int(np.ceil(0.9*msd.shape[0]))
            if not self.inplane and not self.direction :
                ndof = 3
            elif self.inplane:
                ndof = 2
            elif self.direction:
                ndof = 1
            diff = np.polyfit(time[imin:imax],msd[imin:imax],1)[0]/(2.0*ndof) * 10
            print "Diffusion (1e^-5 cm^2/s): %.5f"%diff
        self._write_records()

if __name__ == '__main__' :

    processor = moldyn.TrajectoryProcessor("Calculate the MSD of a group of atoms",
                                            dosubsample=True)
    analysis = MsdAnalysis(processor)
    processor.setup(printargs=True)
    processor.process()
