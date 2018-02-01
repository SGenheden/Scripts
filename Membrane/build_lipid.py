# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to build lipids from a template, similarly to MARTINI INSANE

Is VERY experimental!
"""

import argparse
import os
import xml.etree.ElementTree as ET

import numpy as np

from sgenlib import pdb

class BeadDefinition(object):

    def __init__(self):
        self.name = None
        self.xyz = None

    def parse(self, element):

        if "name" in element.attrib:
            self.name = element.attrib["name"]
        else:
            return

        if "xyz" in element.attrib:
            self.xyz = np.array(element.attrib["xyz"].split(), dtype=float)

    def __str__(self):
        return "%s (%s)"%(self.name,",".join("%.2f"%c for c in self.xyz))

class LipidTemplate(object):

    def __init__(self):
        self.name = None
        self.beads = []
        self.headname = []
        self.tailname = []
        self.head = []
        self.tail = []

    def make(self, bd=3.0):

        struct = pdb.PDBFile()
        res = pdb.Residue()

        for i, bead in enumerate(self.beads):
            atom = pdb.Atom()
            atom.idx = i
            atom.serial = i + 1
            atom.name = bead.name
            atom.resname = self.name
            atom.residue = 1
            atom.set_xyz(bead.xyz*bd)
            res.atoms.append(atom)
            struct.atoms.append(atom)
        struct.residues.append(res)

        allcoord = np.asarray([a.xyz for a in struct.atoms])
        offset = allcoord.mean(axis=0) + 50.0
        for a in struct.atoms:
            a.set_xyz(a.xyz+offset)
        struct.box = np.asarray([100,100,100])

        return struct


    def parse(self, element):

        if "name" in element.attrib:
            self.name = element.attrib["name"]
        else:
            return

        if "head" in element.attrib:
            self.headname = element.attrib["head"].split()
        if "tail" in element.attrib:
            self.tailname = element.attrib["tail"].split()

        for child in element:
            if child.tag != "bead":
                continue
            b = BeadDefinition()
            b.parse(child)
            if b.name is not None:
                self.beads.append(b)
                if b.name in self.headname :
                    self.head.append(b)
                elif b.name in self.tailname :
                    self.tail.append(b)

    def __str__(self):
        return self.name+"\n\t"+"\n\t".join(b.__str__() for b in self.beads)

class LipidCollection(object):

    def __init__(self):
        self.lipids = {}

    def load(self, filename):

        tree = ET.parse(filename)

        # Parse lipids
        for child in tree.getroot():
            if child.tag != "lipid":
                continue
            lipid = LipidTemplate()
            lipid.parse(child)
            if lipid.name is not None:
                self.lipids[lipid.name] = lipid

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Building lipids from templates")
    parser.add_argument('-l','--lipid',help="the lipid to build")
    parser.add_argument('-x','--xml',help="the definition of templates")
    parser.add_argument('-o','--out',help="the output name",default="lipid.pdb")
    parser.add_argument('--bd',type=float,help="the spacing between beads",default=3.0)
    args = parser.parse_args()

    lipidbook = LipidCollection()
    if args.xml is None:
        thispath = os.path.dirname(os.path.abspath(__file__))
        args.xml = os.path.join(thispath,"lipid_templates.xml")
    lipidbook.load(args.xml)

    if args.lipid in lipidbook.lipids:
        struct = lipidbook.lipids[args.lipid].make(bd=args.bd)
        struct.write(args.out)
    else:
        "%s not in the XML file"%args.lipid
