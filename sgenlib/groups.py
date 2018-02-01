# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes and routines to handle atom groups
"""
from __future__ import division, print_function, absolute_import

import sys
import os
import copy
import re
from operator import attrgetter
from ConfigParser import SafeConfigParser

import numpy as np

class AtomGroup :
    """
    Class to store an atom group, the definition
    in terms of atom names and expression in terms
    if atom serial numbers

    Attributes
    ----------
    name : string
      the name of the group
    atoms : list of strings
      the atom names
    serials : dictionary of list of integers
      the serial number for each atom
    expressed : boolean
      if this atom group has been expressed, i.e. if serials has been added
    """
    def __init__(self, name, atoms) :
        self.name = name
        self.atoms = atoms
        self.serials = {a : [] for a in atoms}
        self.expressed = False
    def add_serial(self, atom, serial) :
        """
        Append a serial for a given atom
        """
        try :
            self.serials[atom].append(serial)
            self.expressed = True
        except :
            pass
    def has_atom(self, atom) :
        """
        Check if this group express given atom
        """
        return atom in self.atoms
    def is_consistent(self) :
        """
        Check that all atoms has been expressed equal times
        """
        if len(self.atoms) == 0 : return None

        nadded = len(self.serials[self.atoms[0]])
        for atom in self.atoms[1:] :
            if len(self.serials[atom]) != nadded :
                print(nadded,len(self.serials[atom]),atom)
                return -1
        return nadded
    def write_ndx(self, f, cutoff=20, split=False) :
        """
        Write a Gromac-ndx group
        """
        if not self.expressed : return
        for atom in self.atoms :
            if split : f.write("[%s]\n"%atom)
            for i,serial in enumerate(self.serials[atom],1) :
                f.write("%6d"%serial)
                if i%cutoff == 0 : f.write("\n")
            f.write("\n")
    def __str__(self) :
        return self.name+"="+" ".join(a for a in self.atoms)

class ResidueGroup :
    """
    Class to store a residue group, i.e. a collection
    of AtomGroup objects

    Attributes
    ----------
    name : string
      the residue name
    groups : list of AtomGroup objects
      the atom definitions
    """
    def __init__(self) :
        self.name = ""
        self.groups = []
        self._transmat = None
        self.expanded = False

    def check_consistency(self) :
        """
        Check that all atom groups express equal
        """
        # Find the first non-empy group and store that as a reference
        for off,g in enumerate(self.groups) :
            nadded0 = g.is_consistent()
            if nadded0 is not None : break
        if nadded0 == -1 :
            raise Exception("Unequal expression in %s/%s"%(self.name,self.groups[0].name))

        # Check all other groups
        for g in self.groups[1+off:] :
            nadded = g.is_consistent()
            if nadded == -1 :
                raise Exception("Unequal expression in %s/%s"%(self.name,g.name))
            elif nadded != nadded0 :
                raise Exception("Not all groups in %s are expressed equally"%self.name)

    def indices(self, name):
        """
        Return the indices of the specified name
        """
        return np.asarray([i for i, group in enumerate(self.groups) if group.name == name])

    def expand(self, ncopies, force=False) :
        """
        Make ncopies of the groups
        """
        if self.expanded and not force : return
        newgroups = []
        for i in range(ncopies):
            newgroups.extend([copy.deepcopy(group) for group in self.groups])
        self.groups = newgroups
        self.expanded = True

    def parse(self, section, parser) :
        """
        Parse a ConfigParser section to read in the AtomGroup objects
        """
        self.name = section
        for option in parser.options(section) :
            atoms = parser.get(section,option).split()
            self.groups.append(AtomGroup(option,atoms))

    def transmat(self, atomnames, usechache=True):
        """
        Constructs a transformation matrix from the group names to the group atoms
        """
        if usechache and self._transmat is not None:
            return self._transmat
        self._transmat = np.zeros([len(self.groups), len(atomnames)])
        for i, group in enumerate(self.groups):
            one_over_n = 1 / float(len(group.atoms))
            for atom in group.atoms:
                idx = atomnames.index(atom)
                self._transmat[i, idx] = one_over_n
        return self._transmat

    def write_ndx(self, f) :
        """
        Write Gromacs-ndx group
        """
        for g in self.groups :
            if not g.expressed : continue
            f.write("[%s_%s]\n"%(self.name,g.name))
            g.write_ndx(f)
    def __str__(self) :
        s = "[%s]\n"%self.name
        astr = "\n".join("\t%s"%g for g in self.groups)
        return s+astr

def read_groups(filename) :
    """
    Parse groups from file

    Parameters
    ----------
    filename : string
      the name of the file to parse

    Returns
    -------
    dictionary of ResidueGroup objects
      the parsed groups
    """
    parser = SafeConfigParser()
    parser.read(filename)

    groups = {}
    for section in parser.sections() :
        groups[section] = ResidueGroup()
        groups[section].parse(section,parser)

    return groups

def read_ndxfile(filename) :
    """
    Parse a Gromacs ndx file

    Parameters
    ----------
    filename : string
      the name of the file to parse

    Returns
    -------
    dictionary of lists
      the parsed groups
    """

    ndxstr = ""
    with open(filename, "r") as f :
        ndxstr = "".join(f.readlines())

    ndxsplit = re.split("(\[.+\])", ndxstr)
    ndxdict = {}
    for namestr, indexstr in zip(ndxsplit[1::2], ndxsplit[2::2]):
        ndxdict[namestr[1:-1].strip()] =  \
            map(int, re.sub("\s+"," ",indexstr).strip().split())
    return ndxdict
