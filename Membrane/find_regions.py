# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Script to find the regions of a membrane
"""

import argparse

import numpy as np

from sgenlib import parsing

def _find_intersect(dens1, dens2, z):

    for d1, d2, z in zip(dens1, dens2, z):
        if d1 < d2 :
            return z
    return z[-1]

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Find the regions of a membrane")
    parser.add_argument('-w', '--water', help="the water density", default="w_rdf.txt")
    parser.add_argument('-p', '--phosphate', help="the phosphate density", default="p_rdf.txt")
    parser.add_argument('-g', '--glycerol', help="the glycerol density", default="g_rdf.txt")
    parser.add_argument('-t', '--tail', help="the tail density", default="t_rdf.txt")
    args = parser.parse_args()

    wdens = parsing.parse2ndarray(args.water)
    pdens = parsing.parse2ndarray(args.phosphate)
    gdens = parsing.parse2ndarray(args.glycerol)
    tdens = parsing.parse2ndarray(args.tail)

    n = int(np.round(pdens.shape[0]*0.5))
    leftmaxi = np.argmax(pdens[:n,2])
    rightmaxi = np.argmax(pdens[n:,2])+n
    midz = pdens[leftmaxi,0]+0.5*(pdens[rightmaxi,0]-pdens[leftmaxi,0])
    try :
        midi = np.where(pdens[:,0]==midz)[0][0]
    except:
        midi = np.where(pdens[:,0]==midz)[0]
    if len(midi) == 0 :
        midi = int(leftmaxi + 0.5*(rightmaxi-leftmaxi))
    print "Midz = %.3f"%midz

    wp_left = _find_intersect(wdens[:midi,2], pdens[:midi,2], wdens[:midi,0])
    pg_left = _find_intersect(pdens[:midi,2], gdens[:midi,2], pdens[:midi,0])
    gt_left = _find_intersect(gdens[:midi,2], tdens[:midi,2], gdens[:midi,0])
    gt_right = _find_intersect(gdens[midi:,2][::-1], tdens[midi:,2][::-1], gdens[midi:,0][::-1])
    pg_right = _find_intersect(pdens[midi:,2][::-1], gdens[midi:,2][::-1], pdens[midi:,0][::-1])
    wp_right = _find_intersect(wdens[midi:,2][::-1], pdens[midi:,2][::-1], wdens[midi:,0][::-1])
    print wp_left, pg_left, gt_left, gt_right, pg_right, wp_right
