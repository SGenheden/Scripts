# Author: Samuel Genheden samuel.genheden@gmail.com

import csv
import sys
from collections import namedtuple

import numpy as np

import dblib

ZhangEntry = namedtuple("ZhangEntry",["Solvent","SoluteName","Exper_","katritzky",
                        "katritzky_Difference","COSMO_RS","COSMO_RS_Difference",
                        "TI","TI_error","TI_Difference","Simulation_time"])

class ZhangDb(dblib.SolvDb):

    def read_db(self,filename,type=None,filehandle=None):
        """
        Read in the database from file

        Parameters
        ----------
        filename : string
            the file to read
        type : string, optional
            ignored
        filehandle : string, optional
            ignored
        """
        self._entries = []
        with open(filename,"r") as csvfile:
            reader = csv.DictReader(csvfile,delimiter=";")
            for row in reader:
                keys = row.keys()
                for key in keys:
                    key2 = key.replace(" ","_")
                    key2 = key2.replace(".","_")
                    key2 = key2.replace("-","_")
                    key2 = key2.replace("solvent","Solvent")
                    key2 = key2.replace("solute","SoluteName")
                    if key2 != key:
                        row[key2] = row[key]
                        del row[key]
                self._entries.append(ZhangEntry(**row))
        del self._entries[0]

# Zhang groups
zhang_groups = {}
zhang_groups["carbonyl compound"] = ["2-hexanone","acetone","acetophenone","benzaldehyde","cyclohexanone","formaldehyde"]
zhang_groups["alcohol"] = ["1-butanol","1-octanol","1-pentanol","2-methyl-2-butanol","2-methyl-2-propanol","ethanol","methanol"]
zhang_groups["ether"] = ["anisole","dibutyl_ether","diisopropyl_ether","dimethyl_ether","morpholine","tetrahydrofuran"]
zhang_groups["amine"] = ["diethylamine","morpholine","n-butylamine","triethylamine"]
zhang_groups["alkyl chloride"] = ["12-dichloroethane","1-chlorobutane","chloroform","dichloromethane"]
zhang_groups["carboxylic acid deriv."] = ["ethyl_acetate","methyl _acetate","N-methylformamide","NNdimethylacetamide","NN-dimethylformamide"]
zhang_groups["carboxylic acid ester"] = ["ethyl_acetate","methyl_acetate"]
zhang_groups["carboxylic acid amide"] = ["N-methylformamide","NN-dimethylacetamide","NNdimethylformamide"]
zhang_groups["nitrile"] = ["acetonitrile","benzonitrile"]
zhang_groups["nitro compound"] = ["nitrobenzene","nitromethane"]
zhang_groups["aromatic compound"] = ["2-methylpyridine","3-methylpyridine","acetophenone","anisole","benzaldehyde","benzonitrile","ethylbenzene","isopropylbenzene","nitrobenzene","o-xylene","phenol","pyridine","thiophene","toluene"]
zhang_groups["heterocyclic compound"] = ["2-methylpyridine","3-methylpyridine","morpholine","pyridine","tetrahydrofuran","thiophene"]

db = ZhangDb(sys.argv[1])
"""for entry in db :
    groups = []
    for group,mols in dblib.zhang_groups.iteritems():
        if entry.SoluteName in mols : groups.append(group)
    if not groups : print "***",
    print entry.SoluteName,":"," ".join(groups)

"""
#with open("tmp","r") as f :
#    for line in f.readlines():
#        cols = line.strip().split(", ")
#        print "zhang_groups[\"%s\"] = [%s]"%(cols[0],",".join('\"%s\"'%c.replace(",","") for c in cols[1:]))

print ""

thegroup = "carbonyl compound"
thegroup = "nitro compound"
thegroup = "amine"

for entry in db:
    groups = []
    for group,mols in zhang_groups.iteritems():
        if entry.SoluteName in mols : groups.append(group)
    if groups:
        if thegroup in groups:
            print "%.5f 1"%np.abs(float(entry.TI_Difference))
        else:
            print "%.5f 0"%np.abs(float(entry.TI_Difference))
