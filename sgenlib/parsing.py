# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Routines to parse files and streams
"""

import fileinput
import re

import numpy as np

from units import *

def stdin2ndarray() :
    """
    Parse columns and rows of data from standard input into a numpy.ndarray
    """
    data = []
    for line in fileinput.input() :
        if line[0] not in ["#","@"] :
            data.append(line.strip().split())

    return np.array(data,dtype=float)

def parse2ndarray(filename):
    """
    Parse columns and rows of data from a file into a numpy.ndarray

    Treats lines starting with # and @ as comments
    """
    data = []
    with open(filename,"r") as f :
        data = [line.strip().split() for line in f.readlines()
                 if line[0] not in ["#","@"] ]
    # Check for missing data at the end due to error
    cols0 = len(data[0])
    data = [line for line in data if len(line) == cols0]
    return np.array(data,dtype=float)

def parse_densityfile(filename) :
    """
    Parse densities from g_density
    """
    lines = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while lines[i].startswith("#") : i = i + 1

    legends = []
    while lines[i].startswith('@'):
        if lines[i].find("yaxis  label") > 0 :
            cols = lines[i].split('"')
            ylabel = cols[-2]
        elif re.search(" s\d+ legend",lines[i]):
            cols = lines[i].split('"')
            legends.append(cols[-2])
        i = i + 1

    data = []
    while i < len(lines):
        data.append(lines[i].strip().split())
        i = i + 1
    data = np.array(data, dtype=float)

    densities = {l : col for l, col in zip(legends, data[:,1:].T)}
    return data[:,0], densities, ylabel

##
## Routines to read output from Umbrella Sampling
###

class UmbrellaResults(object) :
    """
    Parent class to encapsulate a results file from an umbrella/zconst
    sampling simulation

    Internally, irrespectively of MD code, the units should
    be Angstroms for lengths and kcal/mol for energy

    Attributes
    ----------
    center : float
        the center of the umbrella potential
    samples : numpy array
        the samples from the simulation
    weight : float
        the weight of the umbrella potential
    """
    def __init__(self,center,weight=None,filename=None,**kwargs) :
        self.samples = None
        self.center = center
        self.weight = weight
        if filename is not None : self.read(filename,**kwargs)
    def __add__(self,other) :
        """ Adds the samples of another simulation onto this one
        """
        both = UmbrellaResults(self.center,self.weight)
        both.samples = np.zeros(self.samples.shape[0]+other.samples.shape[0])
        both.samples[:self.samples.shape[0]] = self.samples
        both.samples[self.samples.shape[0]:] = other.samples
        return both
    def __radd__(self,other) :
        both = UmbrellaResults(other.center,other.weight)
        both.samples = np.zeros(self.samples.shape[0]+other.samples.shape[0])
        both.samples[:self.other.shape[0]] = other.samples
        both.samples[self.other.shape[0]:] = self.samples
        return both
    def read(self,filename,**kwargs) :
        """ Should be implemented by sub-classes
        """
        pass
    def skip(self,nskip,ispart=True) :
        """
        Removes/skip samples at the start of the simulation

        Parameters
        ----------
        nskip : int
            the number or the proportion of samples to skip
        ispart : bool, optional
            if the nskip argument should be interpreted as
            the proportion of the total number of samples
        """
        if self.samples is None or (ispart and nskip < 2): return
        if ispart : nskip = int(np.floor(float(self.samples.shape[0])/float(nskip)))
        self.samples = self.samples[nskip+1:]
    def shorten(self,n,ispart=True) :
        """
        Removes samples at the end of the simulation

        Parameters
        ----------
        n : int
            the number or the proportion of samples to discard
        ispart : bool, optional
            if the n argument should be interpreted as
            the proportion of the total number of samples
        """
        if self.samples is None or (ispart and n < 2): return
        if ispart : n = int(np.floor(float(self.samples.shape[0])/float(n)))
        self.samples = self.samples[:n+1]
    def block_it(self,nblocks) :
        """
        Sub-sample/block the samples

        Parameters
        ----------
        nblocks : int
            the number of blocks

        Returns
        -------
            list of _ResultsFile objects
            the sub-samples/blocks of the original samples
        """
        if nblocks == 1 : return [self]
        blocks = []
        nitems = self.samples.shape[0]
        blocklen = int(np.floor(float(nitems)/float(nblocks)))
        for i in range(nblocks) :
            blocks.append(UmbrellaResults(self.center,self.weight))
            blocks[-1].samples = self.samples[blocklen*i:min(blocklen*(i+1),nitems)]
        return blocks
    def synthesize(self) :
        """
        Replaces the samples with normal distributed data
        according to the mean and standard deviation of the orignal data
        """
        self._samples_orig = self.samples
        mean = self.samples.mean()
        std  = self.samples.std()
        self.samples = random.randn(self.samples.shape[0])*std + mean
    def lower(self) :
        """
        Returns the lower part of the 99% confidence interval
        """
        return self.samples.mean()-2.58*self.samples.std()
    def upper(self) :
        """
        Returns the upper part of the 99% confidence interval
        """
        return self.samples.mean()+2.58*self.samples.std()

class GromacsUmbrellaResults(UmbrellaResults) :
  """
  Read umbrella samples from a Gromacs pullx-files
  """
  def read(self,filename,**kwargs) :
    """
    Read a pullx-file

    Parameters
    ----------
    filename : string
      the name of the file
    **kwargs :
      colidx : the column index to parse
      isforces : if it is a pullf file rather than a pullx file
    """
    data = []
    if "colidx" in kwargs :
      colidx = kwargs["colidx"]
    else :
      colidx = 2
    # The pullf files do not have column 1
    if "isforces" in kwargs and kwargs["isforces"] :
      colidx = colidx - 1

    for line in open(filename,'r').readlines() :
      if not line[0] == "#" and not line[0] == "@" and len(line) > 4 :
        try :
          data.append(line.strip().split()[colidx])
        except :
          pass

    if "isforces" in kwargs and kwargs["isforces"] :
      self.samples =  np.array(data,dtype=float) / 10.0 * ECONV[KJMOL][KCALMOL] # Convert to kcal/mol/A
    else :
      self.samples = np.array(data,dtype=float)*10 # Convert to A
      if colidx != 2 :
        self.samples = -self.samples
      if self.center > 0 : self.samples = np.abs(self.samples)

class PlumedUmbrellaResults(UmbrellaResults) :
  """
  Read umbrella samples from a plumed colvars-file
  """
  def read(self,filename,**kwargs) :
    """
    Read a colvars-file

    Parameters
    ----------
    filename : string
      the name of the file
    **kwargs :
      colidx : the column index to parse
    """

    data = []

    if "colidx" in kwargs :
      colidx = kwargs["colidx"]-1 # To make it compatible with Gromacs umbrella output
    else :
      colidx = 1

    nread = 0
    with open(filename,"r") as f :
      line = f.readline()
      while line :
        if line[0] != "#" :
          data.append(line.strip().split()[colidx])
        line = f.readline()
    self.samples = np.array(data,dtype=float)*10 # Convert to A
    if colidx != 1 :
      self.samples = -self.samples
    if self.center > 0 and "expansion" not in kwargs : self.samples = np.abs(self.samples)

class LammpsUmbrellaResults(UmbrellaResults) :
  """
  Read umbrella samples from a Lammps colvars-file
  """
  def read(self,filename,**kwargs) :
    """
    Read a colvars-file

    Parameters
    ----------
    filename : string
      the name of the file
    **kwargs :
      stride : the stride to parse
    """
    if "stride" in kwargs :
      stride = kwargs["stride"]
    else :
      stride = -1

    data = []
    nread = 0
    with open(filename,"r") as f :
      line = f.readline()
      while line :
        if line[0] != "#" :
          nread = nread + 1
          if stride <= 0 or (stride > 0 and (nread == 1 or (nread-1) % stride == 0)) :
            if len(line.strip().split()) > 1 :
              data.append(line.strip().split()[1])
        line = f.readline()
    self.samples = np.array(data,dtype=float)

class SimpleEnergyFile(UmbrellaResults) :
  """
  Read a simple energy file with total energies
  """
  def read(self,filename,**kwargs) :
    """
    Read a two-columned energy file, where the second column is parsed
    as total energies

    Parameters
    ----------
    filename : string
      the name of the file
    **kwargs :
      ignored at the moment
    """
    data = []
    with open(filename,"r") as f :
      line = f.readline()
      while line :
        if line[0] != "#" :
          data.append(line.strip().split()[1])
        line = f.readline()
    self.samples = np.array(data,dtype=float)
