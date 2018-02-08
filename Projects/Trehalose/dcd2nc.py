# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to convert from CHARMM DCD format to MMTK nc format

Examples
  dcd2nc.py r1_md2.pdb r1_md2_whole.dcd
"""

import sys
import os

from MMTK import *
from nMOLDYN.Core import CHARMMConverter

struct = sys.argv[1]
inname = sys.argv[2]
root, ext = os.path.splitext(os.path.basename(inname))
outname = root + '.nc'

c = CHARMMConverter(struct, inname, outname)
