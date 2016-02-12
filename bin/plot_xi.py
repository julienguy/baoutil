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

parser.add_argument('-d','--data', type = str, default = None, required=True, nargs='*', action='append',
                        help = 'baofit data')
parser.add_argument('-c','--cov', type = str, default = None, required=False,
                        help = 'baofit cov (default is guessed from data)')
parser.add_argument('--mu', type = str, default = None, required=False,
                        help = 'mu range for wedge of the form "mumin:mumax')
parser.add_argument('--rrange', type = str, default = "10:180", required=False,
                        help = 'r range for wedge of the form "rmin:rmax')
parser.add_argument('--rbin', type = float, default = 4.0, required=False,
                        help = 'r bin size')
parser.add_argument('--res', type = str, default = None, required=False, nargs='*', action='append', 
                    help = 'baofit residuals file to plot model')
parser.add_argument('--out', type = str, default = None, required=False,
		                        help = 'output prefix')


args = parser.parse_args()

cov_filename=args.cov



data=[]
for d2 in args.data :
    for d1 in d2 :
       data.append(read_baofit_data(d1))

if cov_filename is None :
    if args.data[0][0].find(".data") :
        cov_filename=string.replace(args.data[0][0],".data",".cov")
    else :
        cov_filename=args.data+".cov"
cov  = read_baofit_cov(cov_filename,n2d=data[0].size,convert=True)

models=[]
if args.res is not None :
    for res2 in args.res :
        for res in res2 :
            mod = read_baofit_model(res,n2d=data[0].size)
            if data[0].size != mod.size :
                print "error data and model don't have same size"
                sys.exit(12)
            models.append(mod)




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

data_colors=["b","r","g","k"]
model_colors=["r","k","k","k"]

for w,wedge in zip(range(nw),wedges) :
    print "plotting mu",wedge
    

    for d,c in zip(data,data_colors) :
        r,xidata,xierr=compute_wedge(d,cov,murange=wedge,rrange=rrange,rbin=args.rbin)        
        scale=r**2
        ax[w].errorbar(r,scale*xidata,scale*xierr,fmt="o",color=c)
    
    for model,c in zip(models,model_colors)  :
        r,ximod,junk=compute_wedge(model,cov,murange=wedge,rrange=rrange,rbin=args.rbin)
        scale=r**2
        ax[w].plot(r,scale*ximod,"-",color=c)
    
    ax[w].set_title(r"$%2.1f < \mu < %2.1f$"%(wedges[w][0],wedges[w][1]))
    ax[w].set_ylabel(r"$r^2 \xi(r)$  [Mpc$^2$ h$^{-2}$]")

plt.xlabel(r"$r$  [Mpc h$^{-1}$]")

plt.show()

if args.out != None:
	plt.savefig(args.out+".png",bbox_inches="tight")

