#!/usr/bin/env python                                                                                                                              
import astropy.io.fits as pyfits
import numpy as np
import sys
import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required=True, help = 'fits file, like exported_cf.fits ')
parser.add_argument('-o','--outfile', type = str, default = None, required=True,help = 'output file')

args = parser.parse_args()

h = pyfits.open(args.infile)
print(h[1].data.dtype.names)


#np.savetxt(args.outfile,xi2d)
#print("wrote",args.outfile)

