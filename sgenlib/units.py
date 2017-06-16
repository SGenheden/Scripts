# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Constants for energy units and conversion, nothing fancy
"""

KJMOL = 1
KCALMOL = 2
KB = {KJMOL : 0.00831446210, KCALMOL : 0.001982923700}
ECONV = {KJMOL : {KJMOL : 1.0, KCALMOL : 1.0/4.184}, KCALMOL : {KJMOL : 4.184, KCALMOL : 1.0}}
ELABEL = {KJMOL : "kJ/mol", KCALMOL : "kcal/mol"}
LABEL2UNIT = {"kJ/mol" : KJMOL, "kcal/mol" : KCALMOL}
