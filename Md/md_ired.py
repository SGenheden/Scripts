# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate S2 order parameters using iRED

By default it calculates S2 of N-H vectors

Examples:
    md_ired.py -f sim.dcd -s ref.pdb
"""

import sys

import numpy as np

from sgenlib import moldyn

class IredAnalysis(moldyn.AnalysisAction):
    """
    Analysis class for iRED

    Attributes
    ----------
    atm1 : MDAnalysis.AtomSelection
        the first atom making up the iRED vector
    atm2 : MDAnalysis.AtomSelection
        the second atom making up the iRED vector
    mat : numpy.ndarray
        the built-up of the correlation matrix
    s2list : list
        the S2 order parameter for each vector at each subsample point
    outname : string
        the name of the output file
    processor : TrajectoryProcessor object
        the trajectory processor calling this analysis
    """
    def __init__(self,processor):
        super(IredAnalysis,self).__init__(processor)
        protsel = processor.universe.selectAtoms(processor.args.pmask)
        self.atm2 = protsel.selectAtoms("name "+processor.args.atoms[1])
        self.atm1 = protsel.selectAtoms("name "+processor.args.atoms[0]+
                                        " and byres name "+processor.args.atoms[1])
        self.mat = np.zeros([len(self.atm1),len(self.atm1)])
        self.s2list = []
        self.outname = processor.args.out
        self.dosubsample = True

    def process(self):
        """
        Building up the correlation matrix, called at each MD snapshot
        """
        v1 = self.atm2.positions-self.atm1.positions
        vlen = 1.0 / np.sqrt((v1*v1).sum(axis=1))

        """mat2 = np.zeros(mat.shape)
        for i in range(nvec):
            for j in range(nvec):
                mat2[j,i] = np.sum(v1[i]*v1[j])*(vlen[i]*vlen[j]) """

        xx1,xx2 = np.meshgrid(v1[:,0],v1[:,0])
        yy1,yy2 = np.meshgrid(v1[:,1],v1[:,1])
        zz1,zz2 = np.meshgrid(v1[:,2],v1[:,2])
        ll1,ll2 = np.meshgrid(vlen,vlen)
        mat0 = (xx1*xx2+yy1*yy2+zz1*zz2)*(ll1*ll2)
        self.mat += 3.0*mat0*mat0-1

    def subsample(self):
        """
        Calculating the S2 order parameters and then zero the correlation matrix
        """
        self.mat = 0.5*(self.mat / float(self.processor.subsamples))
        # Calculating and sorting the eigenvalues and eigenvectors
        eval,evec = np.linalg.eig(self.mat)
        idx = eval.argsort()[::-1]
        eval = eval[idx]
        evec = evec[:,idx]

        prod = evec*evec
        s2 = np.array([1.0-np.sum(eval[5:]*prod[i,5:]) for i in range(prod.shape[0])])

        self.s2list.append(s2)
        self.mat = np.zeros(self.mat.shape)

    def finalize(self):
        """
        Write out the order parameters to disc
        """
        with open(self.outname,"w") as f :
            self.s2list = np.asarray(self.s2list)
            frmstr = "%.5f "*self.s2list.shape[0] # ///frmstr%tuple(rs2)
            for i,(atm,rs2) in enumerate(zip(self.atm1,self.s2list.T)):
                f.write("%d %.5f\n"%(i,rs2.mean()))


if __name__ == '__main__' :

    #print " ".join(sys.argv)
    processor = moldyn.TrajectoryProcessor("Calculate the iRED order parameters",
                    dosubsample=True)
    processor.argparser.add_argument('--atoms',nargs=2,help="the atom names making the vectors",default=["N","H"])
    processor.argparser.add_argument('--pmask',help="the selectiom mask for protein",default="protein")
    processor.argparser.add_argument('-o','--out',help="the output name",default="s2.txt")
    processor.setup()

    analysis = IredAnalysis(processor)
    processor.process()
