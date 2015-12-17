# Author: Samuel Genheden

"""
Cython routines to perform contact analysis on state series

The series are represented by an NxM NumpyArray of uint8
where N is the number of records/points in the series
and M is the number of items
"""

from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.uint8
ctypedef np.uint8_t DTYPE_t
DTYPE_i = np.int
ctypedef np.int_t DTYPE_i_t
DTYPE_f = np.float
ctypedef np.float_t DTYPE_f_t

cdef inline int IMAX(int a,int b): return a if a >= b else b

def lifetime(np.ndarray[DTYPE_t, ndim=2] series) :
    """
    Calculate the life times of the series

    There as M series and the length of the series is N.
    It should only contain 0 and 1, indicating off and on.
    This routine calculate the average and maximum lifetime of the series being on

    Parameters
    ----------
    series : NxM numpy array
        the time series

    Returns
    -------
    M numpy array
        the average lifetime
    M numpy array
        the maximum lifetime
    """
    cdef int nstep = series.shape[0]
    cdef int nitems = series.shape[1]

    cdef np.ndarray[DTYPE_i_t, ndim=1] non      = np.zeros(nitems,dtype=DTYPE_i)
    cdef np.ndarray[DTYPE_i_t, ndim=1] itemsum  = np.zeros(nitems,dtype=DTYPE_i)
    cdef np.ndarray[DTYPE_i_t, ndim=1] itemmax  = np.zeros(nitems,dtype=DTYPE_i)
    cdef np.ndarray[DTYPE_i_t, ndim=1] itemon   = np.zeros(nitems,dtype=DTYPE_i)
    cdef np.ndarray[DTYPE_f_t, ndim=1] itemav   = np.zeros(nitems,dtype=DTYPE_f)

    cdef int i,j

    for i in range(nstep) :
        for j in range(nitems) :
            if series[i,j] == 1 :
                itemon[j] = itemon[j] + 1
            elif  series[i,j] == 0 and itemon[j] > 0 :
                itemsum[j] = itemsum[j] + itemon[j]
                non[j] = non[j] + 1
                itemmax[j] = IMAX(itemmax[j],itemon[j])
                itemon[j] = 0

    for j in range(nitems) :
        if itemon[j] > 0 :
            itemsum[j] = itemsum[j] + itemon[j]
            non[j] = non[j] + 1
            itemmax[j] = IMAX(itemmax[j],itemon[j])
            itemon[j] = 0
        if non[j] > 0 :
            itemav[j] = DTYPE_f(itemsum[j])/DTYPE_f(non[j])

    return itemav,itemmax

def transition(np.ndarray[DTYPE_t, ndim=2] series, np.int nstates) :
    """
    Calculate the transition matrix of the series

    There are M series and the length of the series is N.
    It should only contain integers from 0 to nstates, indicating the series state
    at a particular instance in time.

    Parameters
    ----------
    series : NxM numpy array
        the time series
    nstates : int
        the number of states

    Returns
    -------
    numpy array
        the transition matrix expressed as probabilities
    numpy array
        the raw transition matrix expressed as counts
    """
    cdef int nstep = series.shape[0]
    cdef int nitems = series.shape[1]
    cdef np.ndarray[DTYPE_f_t, ndim=2] M   = np.zeros((nstates+1,nstates+1),dtype=DTYPE_f)
    cdef np.ndarray[DTYPE_f_t, ndim=2] Mcount   = np.zeros((nstates+1,nstates+1),dtype=DTYPE_f)

    cdef int i,j,k,l
    cdef float colsum

    for i in range(nstep-1) :
        for j in range(nitems) :
            k = series[i,j]
            l = series[i+1,j]
            M[k,l] = M[k,l] + 1.0
            Mcount[k,l] = Mcount[k,l] + 1.0

    for i in range(nstates+1) :
        colsum = 0.0
        for j in range(nstates+1) :
            colsum = colsum + M[i,j]
        for j in range(nstates+1) :
            M[i,j] = M[i,j] / colsum

    return M,Mcount

def pairwise_contacts(np.ndarray[DTYPE_t, ndim=2] series) :
    """
    Calculate pairwise contact matrix of the series

    There as M series and the length of the series is N.
    It should only contain 0 and 1, indicating off and on.

    Parameters
    ----------
    series : NxM numpy array
        the time series

    Returns
    -------
    MxM numpy array
        the average pairwise contact number
    """
    cdef int nstep = series.shape[0]
    cdef int nitems = series.shape[1]
    cdef np.ndarray[DTYPE_f_t, ndim=2] M   = np.zeros((nitems,nitems),dtype=DTYPE_f)

    cdef int i,j,k

    for i in range(nstep) :
        for j in range(nitems) :
            if series[i,j] == 1 :
                M[j,j] = M[j,j] + 1
                for k in range(j+1,nitems) :
                    if series[i,k] == 1 :
                       M[j,k] = M[j,k] + 1

    for j in range(nitems) :
        M[j,j] = M[j,j] / float(nstep)
        for k in range(j+1,nitems) :
            M[j,k] = M[j,k] / float(nstep)
            M[k,j] = M[j,k]

    return M

def group_contacts_lead(np.ndarray[DTYPE_t, ndim=2] series, np.ndarray[DTYPE_i_t, ndim=2] groups):
    """
    Calculate join contacts for groups

    The group contacts are determined or-wise. However the leading item needs
    to be on if the whole group contact should be considered as on

    There are M series and the length of the series is N.
    It should only contain 0 and 1, indicating off and on.

    The group array should contain integer indices. Each group is specified
    column-wise, and if the group has less members than the maximum,
    the rest of the indicies should be negative.

    Parameters
    ----------
    series : NxM numpy array
        the time series
    groups : mxG numpy arrays
      the group specifications

    Returns
    -------
    gx1 numpy array
        the  pairwise contact numbers for each group
    """

    cdef int nstep = series.shape[0]
    cdef int nitems = series.shape[1]
    cdef int maxmem = groups.shape[0]
    cdef int ngroups = groups.shape[1]
    cdef np.ndarray[DTYPE_f_t, ndim=1] M   = np.zeros(ngroups,dtype=DTYPE_f)

    cdef int i,g,j
    cdef float on

    for i in range(nstep):
        for g in range(ngroups):
            on = 0.0
            if series[i,groups[0,g]] == 1:
                for j in range(1,maxmem):
                    if groups[j,g] == -1:
                        break
                    if series[i,groups[j,g]] == 1 :
                        on = 1.0
                        break
            M[g] = M[g] + on

    for g in range(ngroups):
        M[g] = M[g] / float(nstep)

    return M

def group_contacts(np.ndarray[DTYPE_t, ndim=2] series, np.ndarray[DTYPE_i_t, ndim=2] groups):
    """
    Calculate join contacts for groups

    The group contacts are determined or-wise.

    There are M series and the length of the series is N.
    It should only contain 0 and 1, indicating off and on.

    The group array should contain integer indices. Each group is specified
    column-wise, and if the group has less members than the maximum,
    the rest of the indicies should be negative.

    Parameters
    ----------
    series : NxM numpy array
        the time series
    groups : mxG numpy arrays
      the group specifications

    Returns
    -------
    gx1 numpy array
        the  pairwise contact numbers for each group
    """

    cdef int nstep = series.shape[0]
    cdef int nitems = series.shape[1]
    cdef int maxmem = groups.shape[0]
    cdef int ngroups = groups.shape[1]
    cdef np.ndarray[DTYPE_f_t, ndim=1] M   = np.zeros(ngroups,dtype=DTYPE_f)

    cdef int i,g,j
    cdef float on

    for i in range(nstep):
        for g in range(ngroups):
            on = 0.0
            for j in range(maxmem):
                if groups[j,g] == -1:
                    break
                if series[i,groups[j,g]] == 1 :
                    on = 1.0
                    break
            M[g] = M[g] + on

    for g in range(ngroups):
        M[g] = M[g] / float(nstep)

    return M

def make_group_series(np.ndarray[DTYPE_t, ndim=2] series, np.ndarray[DTYPE_i_t, ndim=2] groups):

    """
    Make time series for groups

    The group contacts are determined or-wise.

    There are M series and the length of the series is N.
    It should only contain 0 and 1, indicating off and on.

    The group array should contain integer indices. Each group is specified
    column-wise, and if the group has less members than the maximum,
    the rest of the indicies should be negative.

    Parameters
    ----------
    series : NxM numpy array
        the time series
    groups : mxG numpy arrays
      the group specifications


    Returns
    -------
    NxG numpy array
        the time series
    """

    cdef int nstep = series.shape[0]
    cdef int nitems = series.shape[1]
    cdef int maxmem = groups.shape[0]
    cdef int ngroups = groups.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=2] M   = np.zeros((nstep,ngroups),dtype=DTYPE)

    cdef int i,g,j
    cdef int on

    for i in range(nstep):
        for g in range(ngroups):
            on = 0
            for j in range(maxmem):
                if groups[j,g] == -1:
                    break
                if series[i,groups[j,g]] == 1 :
                    on = 1
                    break
            M[i,g] = on

    return M

def make_group_series_lead(np.ndarray[DTYPE_t, ndim=2] series, np.ndarray[DTYPE_i_t, ndim=2] groups):

    """
    Make time series for groups

    The group contacts are determined or-wise. However the leading item needs
    to be on if the whole group contact should be considered as on.

    There are M series and the length of the series is N.
    It should only contain 0 and 1, indicating off and on.

    The group array should contain integer indices. Each group is specified
    column-wise, and if the group has less members than the maximum,
    the rest of the indicies should be negative.

    Parameters
    ----------
    series : NxM numpy array
        the time series
    groups : mxG numpy arrays
      the group specifications


    Returns
    -------
    NxG numpy array
        the time series
    """

    cdef int nstep = series.shape[0]
    cdef int nitems = series.shape[1]
    cdef int maxmem = groups.shape[0]
    cdef int ngroups = groups.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=2] M   = np.zeros((nstep,ngroups),dtype=DTYPE)

    cdef int i,g,j
    cdef int on

    for i in range(nstep):
        for g in range(ngroups):
            on = 0
            if series[i,groups[0,g]] == 1:
                for j in range(1,maxmem):
                    if groups[j,g] == -1:
                        break
                    if series[i,groups[j,g]] == 1 :
                        on = 1
                        break
            M[i,g] = on

    return M
