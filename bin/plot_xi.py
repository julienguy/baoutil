#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import string

from baoutil.io import read_baofit_data
from baoutil.io import read_baofit_cov
from baoutil.io import read_baofit_model
from baoutil.wedge import compute_wedge

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-d','--data', type = str, default = None, required=True,
                        help = 'baofit data')
parser.add_argument('-c','--cov', type = str, default = None, required=False,
                        help = 'baofit cov (default is guessed from data)')
parser.add_argument('--mu', type = str, default = None, required=False,
                        help = 'mu range for wedge of the form "mumin:mumax')
parser.add_argument('--rrange', type = str, default = "10:180", required=False,
                        help = 'r range for wedge of the form "rmin:rmax')
parser.add_argument('--rbin', type = float, default = 4.0, required=False,
                        help = 'r bin size')
parser.add_argument('--res', type = str, default = None, required=False,
                        help = 'baofit residuals file to plot model')
parser.add_argument('--out', type = str, default = None, required=False,
		                        help = 'output prefix')


args = parser.parse_args()

cov_filename=args.cov

if cov_filename is None :
    if args.data.find(".data") :
        cov_filename=string.replace(args.data,".data",".cov")
    else :
        cov_filename=args.data+".cov"
data = read_baofit_data(args.data)        
cov  = read_baofit_cov(cov_filename,n2d=data.size,convert=True)

if args.res is not None :
    mod = read_baofit_model(args.res,n2d=data.size)
    if data.size != mod.size :
        print "error data and model don't have same size"
        sys.exit(12)
else :
    mod = None



if args.mu :
    try :
        vals=string.split(args.mu,":")
        if len(vals)!=2 :
            print "incorrect format for mu range '%s', expect mumin:mumax"%args.mu
            sys.exit(12)
        mumin=string.atof(vals[0])
        mumax=string.atof(vals[1])
        wedges=[[mumin,mumax]]
    except ValueError,e:
        print e
        print "incorrect format for mu range '%s', expect mumin:mumax"%args.mu
        sys.exit(12)
else :
    wedges= [[0.8,1.0],[0.5,0.8],[0.0,0.5]]

if args.rrange :
    try :
        vals=string.split(args.rrange,":")
        if len(vals)!=2 :
            print "incorrect format for r range '%s', expect rmin:rmax"%args.rrange
            sys.exit(12)
        rmin=string.atof(vals[0])
        rmax=string.atof(vals[1])
        rrange=[rmin,rmax]
    except ValueError,e:
        print e
        print "incorrect format for r range '%s', expect rmin:rmax"%args.rrange
        sys.exit(12)
else :
    rrange=[10,180]


nw=len(wedges)
if nw > 1 :
    f, ax = plt.subplots(nw, sharex=True, sharey=False)
else :
    f=plt.figure()
    ax=[]
    ax.append(plt.subplot(1,1,1))



for w,wedge in zip(range(nw),wedges) :
    print "plotting mu",wedge
    
    r,xidata,xierr=compute_wedge(data,cov,murange=wedge,rrange=rrange,rbin=args.rbin)
    scale=r**2
    ax[w].errorbar(r,scale*xidata,scale*xierr,fmt="o",color="b")

    if mod is not None :
        r2,ximod,junk=compute_wedge(mod,cov,murange=wedge,rrange=rrange,rbin=args.rbin)
        ax[w].plot(r,scale*ximod,"-",color="r")
    
    ax[w].set_title(r"$%2.1f < \mu < %2.1f$"%(wedges[w][0],wedges[w][1]))
    ax[w].set_ylabel(r"$r^2 \xi(r)$  [Mpc$^2$ h$^{-2}$]")

plt.xlabel(r"$r$  [Mpc h$^{-1}$]")

plt.show()

if args.out != None:
	plt.savefig(args.out+".png",bbox_inches="tight")

