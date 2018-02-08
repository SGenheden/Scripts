# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Module to read in the Minnesota Solvation database
and then supply generators that can be used to iterate over
the entries
"""

import csv
from collections import namedtuple
import re


Entry = namedtuple("Entry",['No', 'FileHandle', 'SoluteName', 'Formula', 'Subset',
                            'Charge', 'Level1', 'Level2', 'Level3', 'Solvent',
                            'DeltaGsolv', 'type', 'eps', 'n', 'alpha', 'beta',
                            'gamma', 'phi2', 'psi2', 'beta2', 'H', 'C', 'HC',
                            'CC', 'CC2', 'N', 'HN', 'HN2', 'CN', 'NC', 'NC2', 'NC3',
                            'O', 'HO', 'HO2', 'OC', 'CO2', 'ON', 'OO', 'F', 'Cl', 'Br',
                             'I', 'FC', 'ClC', 'BrC', 'IC', 'Si', 'OSi', 'P', 'HP', 'OP',
                            'S', 'HS', 'OS', 'SP', 'SS', 'TotalArea'])

class SolvDb(object):
    """
    Class to read the database and to provide generators
    over the entries. The entries themselves are hidden
    """
    def __init__(self,filename=None,type=None,filehandle=None):
        self._entries = []
        if filename is not None : self.read_db(filename,type,filehandle)

    def __len__(self):
        return len(self._entries)

    def __getitem__(self,key):
        if key > len(self._entries) :
            raise IndexError("Outside number of entries")
        return self._entries[key]

    def itersolvent(self,solvent):
        """
        Generator that gives all entries for a particular solvent
        """
        for entry in self._entries:
            if entry.Solvent == solvent : yield entry

    def itersolutelist(self,solvent,solutes):
        """
        Generator that gives all entries for a particular solvent and
        where the solute is in a list of solutes
        """
        for entry in self._entries:
            if entry.Solvent == solvent and entry.SoluteName in solutes :
                yield entry

    def itersoluteoverlap(self,solvent1,solvent2):
        """
        Generator that gives all entries where the solute is found in two
        different solvents
        """
        solutes2 = [entry for entry in self.itersolvent(solvent2)]
        solutelist = [entry.SoluteName for entry in solutes2]
        for entry1 in self.itersolutelist(solvent1,solutelist):
            yield entry1,solutes2[solutelist.index(entry1.SoluteName)]

    def read_db(self,filename,type=None,filehandle=None):
        """
        Read in the database from file

        Parameters
        ----------
        filename : string
            the file to read
        type : string, optional
            if given, this restricts reading to entries with this type
        filehandle : string, optional
            if given, this restricts reading to entries that match this regular expression
        """
        self._entries = []
        with open(filename,"r") as csvfile:
            reader = csv.DictReader(csvfile,delimiter="\t")
            for row in reader:
                if type is not None and row["type"] != type : continue
                if filehandle is not None and not re.match(filehandle,row["FileHandle"]) : continue
                # This ugly piece of code removes non-allowed characters from the keys
                for lbl in [("No","No."),("phi2","phi**2"),("psi2","psi**2"),("beta2","beta**2")]:
                    row[lbl[0]]=row[lbl[1]]
                    del row[lbl[1]]
                self._entries.append(Entry(**row))

supergroups = ['alcohol',
                'aldehyde',
                'alkane',
                'alkene',
                'alkyl bromide',
                'alkyl chloride',
                'amine',
                'aromatic compound',
                'carbonitrile',
                'carboxylic acid',
                'carboxylic acid ester',
                'ether',
                'halogen derivative',
                'heterocyclic compound',
                'ketone',
                'nitro compound',
                'phenol or hydroxyhetarene']

def transform_groups(groups,entry):
    """Transform a set of groups to a list of supergroups"""
    groups2 = []
    for group in groups:
        group2 = group.replace("primary","")
        group2 = group2.replace("secondary","")
        group2 = group2.replace("tertiary","")
        if group2.find("ether") > -1:
            group2 = "ether"
        elif group2.find("aliphatic amine") > -1:
            group2 = "amine"
        group2 = group2.strip()
        if group2 in supergroups:
            groups2.append(group2)
    if not groups2 and entry.SoluteName.find("cyclo") == -1:
        groups2.append("alkane")
    return groups2
