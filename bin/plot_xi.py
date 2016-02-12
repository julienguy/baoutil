#!/usr/bin/env python
import numpy as np
import pylab
import sys
import argparse
import string

from baoutil.io import read_baofit_data
from baoutil.io import read_baofit_cov
from baoutil.io import read_baofit_model
from baoutil.wedge import compute_wedge

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--data', type = str, default = None, required=True,
                        help = 'baofit data')
parser.add_argument('--cov', type = str, default = None, required=False,
                        help = 'baofit cov (default is guessed from data)')
parser.add_argument('--res', type = str, default = None, required=False,
                        help = 'baofit residuals file to plot model')

args = parser.parse_args()

cov_filename=args.cov

if cov_filename is None :
    if args.data.find(".data") :
        cov_filename=string.replace(args.data,".data",".cov")
    else :
        cov_filename=args.data+".cov"
data = read_baofit_data(args.data)        
cov  = read_baofit_cov(cov_filename)
if args.res is not None :
    mod = read_baofit_model(args.res,n2d=data.size)
    if data.size != mod.size :
        print "error data and model don't have same size"
        sys.exit(12)
else :
    mod = None

rrange=[10,180]
wedges= [[0.8,1.0],[0.5,0.8],[0.0,0.5]]
rbin=4

nw=len(wedges)
if nw > 1 :
    f, ax = pylab.subplots(nw, sharex=True, sharey=False)
else :
    f=pylab.figure()
    ax=[]
    ax.append(pylab.subplot(1,1,1))



for w,wedge in zip(range(nw),wedges) :
    print "plotting mu",wedge
    
    r,xidata,xierr=compute_wedge(data,cov,murange=wedge,rrange=rrange,rbin=rbin)
    scale=r**2
    ax[w].errorbar(r,scale*xidata,scale*xierr,fmt="o",color="b")

    if mod is not None :
        r2,ximod,junk=compute_wedge(mod,cov,murange=wedge,rrange=rrange,rbin=rbin)
        ax[w].plot(r,scale*ximod,"-",color="r")
    
    ax[w].set_title(r"$%2.1f < \mu < %2.1f$"%(wedges[w][0],wedges[w][1]))
    ax[w].set_ylabel(r"$r^2 \xi(r)$  [Mpc$^2$ h$^{-2}$]")

pylab.xlabel(r"$r$  [Mpc h$^{-1}$]")

pylab.show()

