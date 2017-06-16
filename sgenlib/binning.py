# Author: Samuel Genheden samuel.genheden@gmail.com
"""
Routines for binning data, histogramming and such
"""

import numpy as np

def make_bins(data,nbins,boundaries=None) :
    """
    Make bin edges for histogramming

    Parameters
    ----------
    data : numpy array
        the data to be histogrammed
    nbins : int
        the number of bins to create
    boundaries : list of float, optional
        the minimum and maximum of the edges,
        if not supplied used the 99% confidence interval of the data (assuming normality)

    Results
    -------
    numpy array
        the bins edges
    """
    def get_99_conf(data) :
        mean = data.mean()
        std  = data.std()
        return (mean-2.58*std,mean+2.58*std)

    if boundaries == None :
        mn,ma = get_99_conf(data[0])
        for d in  data[1:] :
            lo,hi = get_99_conf(d)
            mn = min(mn,lo)
            ma = max(ma,hi)
    else :
        mn, ma = boundaries

    return np.linspace(mn, ma, nbins+1, endpoint=True)


def make_histograms(data,bins,boundaries=None) :
    """
    Histogram set of data using the same edges

    Mainly used by UmbrellaSim class, but could be useful in other circumstances

    Parameters
    ----------
    data : numpy array
        the data to be histogrammed
    bins : int or numpy array
        the number of bins to use or the edges of the histogram
    boundaries : list of float, optional
        the minimum and maximum of the edges

    Returns
    -------
    list of numpy array
        the histogram of the data
    numpy array
        the bin edges
    """
    if not np.iterable(bins) :
        bins = make_bins(data,bins,boundaries)

    histograms = [np.histogram(d,bins)[0] for d in data]
    return histograms,bins

def make_reweighted_histograms(data,weights,bins,boundaries=None) :
    """
    Weight-histogram set of data using the same edges

    Mainly used by UmbrellaSim class, but could be useful in other circumstances

    Parameters
    ----------
    data : numpy array
        the data to be histogrammed
    weights : numpy array
        the weight of the data
    bins : int or numpy array
        the number of bins to use or the edges of the histogram
    boundaries : list of float, optional
        the minimum and maximum of the edges

    Returns
    -------
    list of numpy array
        the histogram of the data
    numpy array
        the bin edges
    """
    if not np.iterable(bins) :
        bins = make_bins(data,bins,boundaries)

    histograms = [np.histogram(d,bins,weights=w)[0] for d,w in zip(data,weights)]
    return histograms,bins
