
import sys
import os

import parmed

ambname = sys.argv[1]
crdname = os.path.splitext(ambname)[0]+".crd"
mol2name = os.path.splitext(ambname)[0]+".mol2"
pdbname = os.path.splitext(ambname)[0]+".pdb"

top = parmed.amber.AmberParm(ambname, xyz=crdname)
top.save(mol2name, format="mol2")
top.save(pdbname, format="pdb")
