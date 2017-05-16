# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Numerical algorithms
"""
import numpy as np

def trapz(x, av, stds) :
    """
    Trapezoid integration with error propagation
    """
    weights = np.zeros(x.shape[0])
    weights[0] = 0.5*(x[0]+x[1])
    weights[-1] = 1.0 - 0.5*(x[-1]+x[-2])
    weights[1:-1] = 0.5 * (x[2:]-x[:-2])
    var = (weights**2*stds**2).sum()
    return np.trapz(av, x), np.sqrt(var)
