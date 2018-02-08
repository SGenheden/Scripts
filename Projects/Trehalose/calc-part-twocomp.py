#!/Users/samuel/anaconda/envs/oldnumpyenv/bin/python
# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse

from MMTK import *
from nMOLDYN.Core.Logger import LogCreate
from nMOLDYN.Analysis.Templates import StaticCoherentStructureFactor_serial

LogCreate(['console', 'file'])

argparser = argparse.ArgumentParser(description="Script calculate structure factors for two-component system")
argparser.add_argument('-f', '--file', help="the input trajectory")
argparser.add_argument('-o', '--out',  help="the output file")
argparser.add_argument('-d','--deuteration',choices=['none','water','trehalose','trehalosex', 'bothx'],help="the type of deuteration", default="none")
args = argparser.parse_args()

parameters = {}

parameters['version'] = "3.0.10"
parameters['estimate'] = "no"
parameters['pyroserver'] = "monoprocessor"
parameters['qunits'] = "nm^-1"
parameters['qvectors'] = {'qgeometry': 'spatial', 'qshellwidth': 1.0, 'qshellvalues': '0.0:100.0:1.0', 'hkls': None, 'qvectorspershell': 50}
parameters['subset'] = "object objectname * atomelement hydrogen"
parameters['timeinfo'] = "1:20001:5"
parameters['weights'] = "coherent"

if args.deuteration == "none":
    parameters['deuteration'] = None
elif args.deuteration == "water":
    parameters['deuteration'] = "object objectname SOL atomelement hydrogen"
elif args.deuteration == "trehalose":
    parameters['deuteration'] = "object objectname 0GA atomelement hydrogen OR objectname 1GA atomelement hydrogen"
elif args.deuteration == "trehalosex":
    parameters['deuteration'] = "object objectname 0GA atomname 0GA.H2O,0GA.H3O,0GA.H4O,0GA.H6O OR objectname 1GA atomname 1GA.H2O,1GA.H3O,1GA.H4O,1GA.H6O"
elif args.deuteration == "bothx":
    parameters['deuteration'] = "object objectname 0GA atomname 0GA.H2O,0GA.H3O,0GA.H4O,0GA.H6O OR objectname 1GA atomname 1GA.H2O,1GA.H3O,1GA.H4O,1GA.H6O or objectname SOL atomelement hydrogen"

parameters['trajectory'] = args.file
parameters['output'] = args.out

analysis = StaticCoherentStructureFactor_serial(parameters)
analysis.runAnalysis()
