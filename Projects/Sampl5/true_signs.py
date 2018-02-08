# Author: Samuel Genheden samuel.genheden@gmail.com

import sys

import numpy as np

from sgenlib import parsing

if __name__ == '__main__' :

    data = parsing.parse2ndarray(sys.argv[1])
    true_signs = (np.abs(data[:,0])>=1.96*data[:,1]).sum()
    print "Number of true signs is %d of %d, that is %d percent"%(true_signs,
        data.shape[0], float(true_signs)/float(data.shape[0])*100)
    print "Predictions with false signs %s"% \
        ", ".join("%.1f+-%.2f"%tuple(d) for d in data[np.abs(data[:,0])<1.96*data[:,1],:])
