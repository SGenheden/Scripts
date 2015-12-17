# Author: Samuel Genheden, samuel.genheden@gmail.com

"""
Program to output a netcdf xy variables
"""

import argparse

from scipy.io import netcdf


if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description="Output NetCDF variables")
    parser.add_argument('file',help="the file")
    parser.add_argument('-x','--x',help="the x variable")
    parser.add_argument('-y','--y',help="the y variable")
    args = parser.parse_args()

    with netcdf.netcdf_file(args.file, 'r') as f:
        #print f.variables.keys()
        for x,y in zip(f.variables[args.x],f.variables[args.y]):
            print x,y
