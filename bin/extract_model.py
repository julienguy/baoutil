#!/usr/bin/env python                                                                                                                              
import h5py
import numpy as np
import sys
import argparse
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required=True, help = 'h5 file, like fit_lya_lyaxlya_lya.h5')
parser.add_argument('-o','--outfile', type = str, default = None, required=True,help = 'output file')
parser.add_argument('-w','--what', type = str, default = None, required=False,help = 'what, like "LYA(LYA)-LYA(LYA)"')

args = parser.parse_args()

f = h5py.File(args.infile)
names = [k for k in f]
print("names are",names)
if args.what is None : 
    args.what = names[0]
    print("will show:",args.what)



for i,k in enumerate(f):                                                                            
    if k != args.what: continue                                                           
    xi2d = np.array(f[str(k)][u'fit'])
    print(xi2d.shape)
    np.savetxt(args.outfile,xi2d)
    print("wrote",args.outfile)

