# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to help with the processing of MD trajectories
"""

import argparse
import sys
from collections import namedtuple

import MDAnalysis as md

class TrajectoryProcessor(object):
    """
    Class to process an MD trajectory

    The program that uses this initialises an instance that setups up an
    argparse command-line interpreter.
    The program can add its own arguments to this parser before calling setup()
    that parse the arguments and initialises the md universe object

    To analyse the trajectory the program then calls the process() routine that
    calles AnalysisAction objects in the analyses list. These objects needs to
    be initiated before a call to process().

    Attributes
    ----------
    analyses : list of AnalysisAction objects
        the analyses to call
    argparser : argparse.ArgumentParser
        the command line argument parser
    args : argparse.Namespace
        the parsed arguments
    currtime : float
        the current time in ps when processing the trajectory
    currbox : numpy.ndarray
        the box of the current snapshot
    currsnap : MDAnalysis.TimeStep object
        the current snapshot
    dosubsample : boolean
        if analyses doing subsampling
    dt : float
        the number of ps per snapshot
    every : int
        the process frequency
    freq : int
        the subsample frequency
    skip : float
        the number of ps to skip
    subsamples : int
        the number of subsamples
    universe : MDAnalysis.Universe
        the MDAnalysis universe holding the trajectory and its state
    """


    def __init__(self, description,dosubsample=False):
        self.argparser = argparse.ArgumentParser(description=description)
        self.argparser.add_argument('-f', '--file', nargs="+", help="the trajectory file.")
        self.argparser.add_argument('-s', '--struct', help="a structure file")
        self.argparser.add_argument('--skip', type=int, help="skip this many snapshots", default=0)
        self.argparser.add_argument('--every',type=int,help="print out every X frame",default=1)
        self.argparser.add_argument('--dt', type=float, help="the number of ps for each snapshot", default=10)
        if dosubsample:
            self.argparser.add_argument('--freq',type=int,help="the analysis frequency",default=100000)
        self.universe = None
        self.dosubsample=dosubsample
        self.currtime = 0
        self.dt = 0
        self.skip = 0
        self.every = 0
        self.subsamples = 0
        self.freq = 0
        self.analyses = []

    def setup(self):
        """
        Parsing command-line arguments and setting up the MD universe
        This routine should be called before process()
        """
        self.args = self.argparser.parse_args()
        self.skip = self.args.skip
        self.dt = self.args.dt
        self.every = self.args.every
        if self.dosubsample:
            self.freq = self.args.freq
            self.subsamples = self.freq/self.dt
        self.setup_universe()

    def setup_universe(self):
        if self.args.file is None or self.args.struct is None :
            raise Exception("Both a universe and structure file needs to be specified")
        self.universe = md.Universe(self.args.struct, self.args.file)

    def process(self):
        """
        Loop over the MD trajectory and process it, calling the objects in
        the analyses list
        """
        if self.universe is None or not self.analyses: return

        for ti, ts in enumerate(self.universe.trajectory,1):
            self.currsnap = ts
            self.currbox = ts.dimensions[:3]
            self.currtime = ti*self.dt
            if self.currtime % 10000 == 0:
                print "%d" % (self.currtime),
                sys.stdout.flush()
            if self.currtime < self.skip:
                continue
            if ti % self.every != 0 :
                continue
            for analysis in self.analyses:
                analysis.process()
                if self.dosubsample and analysis.dosubsample and \
                        self.currtime % self.freq == 0:
                    analysis.subsample()

        print ""
        for analysis in self.analyses :
            analysis.finalize()

MDRecord = namedtuple("MDRecord",["time","value"])

class AnalysisAction(object):
    """
    Base class for analysis actions that are called by a TrajectoryProcessor
    object. Classes that wants to implement analysis a particular action
    should inherit from this.
    """
    def __init__(self,processor):
        self.processor = processor
        self.dosubsample = False
        processor.analyses.append(self)

    def process(self):
        pass

    def subsample(self):
        pass

    def finalize(self):
        pass

    def _write_records(self,postfix=''):
        """
        Helper routine to write out a list of MDRecord to disc
        """
        if self.records :
            with open(self.processor.args.out+postfix,'w') as f :
                for entry in self.records:
                    f.write("%.0f %.3f\n"%(entry.time,entry.value))
