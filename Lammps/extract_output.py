# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to extract the thermodynamic output from a Lammps output file

Example:
    extract_output.py out.sim > out.sim_data
"""

import sys

with open(sys.argv[1],"r") as f :
    line = f.readline()
    parsed_lines = []
    parsing = False
    while line:
        if line.startswith("Setting up "):
            parsed_lines = [] # This will ensure we are only keeping the last run
            parsing = True
            line = f.readline() # Unit style
            line = f.readline() # Current step
            line = f.readline() # Time step
            line = f.readline() # Memory
            line = f.readline() # Column header line
            parsed_lines.append("# "+line)
        elif line.startswith("Loop time of"):
            parsing = False
        elif parsing:
            parsed_lines.append(line)
        line = f.readline()

for line in parsed_lines:
    print line,
