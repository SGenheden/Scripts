# Author: Samuel Genheden samuel.genheden@gmail.com
"""
This script uses the encore library
"""
import copy
import argparse

import numpy as np
import MDAnalysis as md
import encore
from encore.similarity import harmonic_ensemble_similarity, bootstrap_coordinates

class MyEnsemble(encore.Ensemble) :

    def get_coordinates(self, subset_selection_string=None, firsthalf=True):

        if not subset_selection_string:
            subset_selection_string = self.atom_selection_string
        subset_selection = self.universe.select_atoms(subset_selection_string)

        if len(subset_selection) == 0:
            logging.error("ERROR: selection \'%s\' not found in topology."% subset_selection_string)
            exit(1)
        try:
            subset_coordinates = self.universe.trajectory.timeseries(subset_selection, skip=self.frame_interval, format='fac')
        except:
            n_coordinates = 0
            k = 0
            for i,time_step in enumerate(self.universe.trajectory):
                if (i % self.frame_interval) == 0:
                    n_coordinates += 1
            subset_coordinates = numpy.zeros(tuple([n_coordinates]) + subset_selection.coordinates().shape)

            for i, time_step in enumerate(self.universe.trajectory):
                if (i % self.frame_interval) == 0:
                    subset_coordinates[k] = subset_selection.coordinates(time_step)
                    k+=1

        n = int(0.5*subset_coordinates.shape[0])
        if firsthalf :
            return subset_coordinates[:n,:]
        else :
            return subset_coordinates[n:,:]

if __name__ == '__main__':

    argparser = argparse.ArgumentParser(description="Program to calculate HES metric on the two halves of a trajectory")
    argparser.add_argument('-s','--struct',help="the filename of a PDB file")
    argparser.add_argument('-f','--file',help="the DCD trajectory")
    argparser.add_argument('-o','--out',help="output prefix")
    args = argparser.parse_args()

    nboots = 500

    uni = md.Universe(args.struct, args.file)
    ensemble = encore.Ensemble(uni, trajectory=args.file)
    n = int(0.5*ensemble.coordinates.shape[0])

    ensemble1 = copy.copy(ensemble)
    ensemble1.coordinates = ensemble.coordinates[:n,:]
    bootcoord1 = bootstrap_coordinates(ensemble1.coordinates, nboots)
    sigma1 = encore.covariance_matrix(ensemble1)

    ensemble2 = copy.copy(ensemble)
    ensemble2.coordinates = ensemble.coordinates[n:,:]
    bootcoord2 = bootstrap_coordinates(ensemble2.coordinates, nboots)
    sigma2 = encore.covariance_matrix(ensemble2)

    boots = []
    for c1, c2 in zip(bootcoord1, bootcoord2):
        x1 = np.average(c1, axis=0).flatten()
        x2 = np.average(c2, axis=0).flatten()
        boots.append(harmonic_ensemble_similarity(x1=x1, x2=x2,
                    sigma1=sigma1, sigma2=sigma2))
    boots = np.asarray(boots)
    print "%.3f\t%.3f"%(boots.mean(), boots.std())
