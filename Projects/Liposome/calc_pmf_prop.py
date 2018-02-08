# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Script to calculate PMF properties, i.e. penetration and water/lipid barriers
"""

import sys

import numpy as np

from sgenlib import parsing

data = parsing.parse2ndarray(sys.argv[1])
x = data[:,0]
pmf = data[:,1]
err = data[:,2]

l = data.shape[0]
imid = int(np.floor(0.5*l))
ifourth = int(np.floor(0.25*l))
dt = np.ceil(1/(x[1]-x[0]))
margin = 10

# This is the index of the maximum in the middle
midmaxi = np.argmax(pmf[imid-dt*margin:imid+dt*margin])+imid-dt*margin

# Now we get the index of the minimum in the outer leaflet
outmid = midmaxi + ifourth
outmini = np.argmin(pmf[outmid-dt*margin:outmid+dt*margin])+outmid-dt*margin

# Now we get the index of the minimum in the inner leaflet
inmid = midmaxi - ifourth
inmini = np.argmin(pmf[inmid-dt*margin:inmid+dt*margin])+inmid-dt*margin

dGdepth_out = pmf[outmini] - pmf[-1]
dGpen_out   = pmf[midmaxi] - pmf[outmini]
dGpen_in    = pmf[inmini]  - pmf[midmaxi]
dGdepth_in  = pmf[0]       - pmf[inmini]

errdepth_out = np.sqrt(err[outmini]**2 + err[-1]**2)
errpen_out   = np.sqrt(err[midmaxi]**2 + err[outmini]**2)
errpen_in    = np.sqrt(err[inmini]**2  + err[midmaxi]**2)
errdepth_in  = np.sqrt(err[0]**2       + err[inmini]**2)

print "\t".join("%.3f"%g for g in (dGdepth_out, errdepth_out,
                                    dGpen_out, errpen_out,
                                    dGpen_in, errpen_in,
                                    dGdepth_in, errdepth_in))
