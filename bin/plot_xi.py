#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import string
import os.path

from baoutil.io import read_baofit_data
from baoutil.io import read_baofit_cov
from baoutil.io import read_baofit_fits
from baoutil.io import read_baofit_model
from baoutil.wedge import compute_wedge
from baoutil.wedge import compute_wedge_with_ivar

plt.rcParams["font.family"]="serif"
plt.rcParams["font.size"]=16.0


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-d','--data', type = str, default = None, required=True, nargs='*',
                        help = 'baofit data')
parser.add_argument('--mu', type = str, default = None, required=False,
                        help = 'mu range for wedge of the form "mumin:mumax,mumin:mumax ...')
parser.add_argument('--rrange', type = str, default = "10:180", required=False,
                        help = 'r range for wedge of the form "rmin:rmax')
parser.add_argument('--rbin', type = float, default = 4.0, required=False,
                        help = 'r bin size')
parser.add_argument('--rpmin', type = float, default = 0, required=False,
                        help = 'min r_parallel')
parser.add_argument('--res', type = str, default = None, required=False, nargs='*',
                    help = 'baofit residuals file to plot model')
parser.add_argument('--out', type = str, default = None, required=False,
		                        help = 'output prefix')
parser.add_argument('--chi2', action="store_true",
		help = 'compute chi2 of wedges (if only one data set and model)')

parser.add_argument('--noshow', action="store_true",
		            help = 'prevent the figure window from displaying')
parser.add_argument('--rpower', type = int, default = 2, required=False,
                        help = 'r power for display')
parser.add_argument('--flip', action="store_true",
		            help = 'flip plot (useful for Lya-QSO cross-corr)')

#parser.add_argument('--ivar_weight', action="store_true",
# help = 'use inverse variance to combine the bins')
parser.add_argument('--no_ivar_weight', action="store_true",
                    help = 'do not use inverse variance to combine the bins')


args = parser.parse_args()


data=[]
cov=[]
for d1 in args.data :
        if d1.find(".fits")>=0 :
            d,c=read_baofit_fits(d1)
            data.append(d)
            cov.append(c)
        else :
            d=read_baofit_data(d1)
            data.append(d)
            if d1.find(".data")>=0 :
                cov_filename=string.replace(d1,".data",".cov")
                if not os.path.isfile(cov_filename) :
                        cov_filename=string.replace(d1,".data","-cov.fits")
                if os.path.isfile(cov_filename) :
                        c  = read_baofit_cov(cov_filename,n2d=d.size,convert=True)
                        cov.append(c)
                else :
                    print "warning, cannot find covariance of ",d1
                    cov.append(np.eye(d.size)*1e-12)    
            else :
                print "warning, cannot guess covariance of ",d1
                cov.append(np.eye(d.size)*1e-12)    


models=[]
if args.res is not None :
    for i,res in enumerate(args.res) :
        if i<len(data) :
            d=data[i]
        else :
            d=data[0]
        mod = read_baofit_model(res,n2d=d.size)
        if d.size != mod.size :
            print "error data and model don't have same size"
            sys.exit(12)
        models.append(mod)
    



if args.mu :
    wedges=[]
    try :
        wedge_strings=string.split(args.mu,",")
        for wedge_string in wedge_strings :
            print wedge_string
            vals=string.split(wedge_string,":")
            if len(vals)!=2 :
                print "incorrect format for mu range '%s', expect mumin:mumax"%args.mu
                sys.exit(12)
            mumin=string.atof(vals[0])
            mumax=string.atof(vals[1])
            wedges.append([mumin,mumax])
    except ValueError,e:
        print e
        print "incorrect format for mu range '%s', expect mumin:mumax"%args.mu
        sys.exit(12)
else :
    wedges= [[0.8,1.0],[0.5,0.8],[0.0,0.5]]
    wedges= [[0.95,1.0],[0.8,0.95],[0.5,0.8],[0.0,0.5]]

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

if args.flip :
    ax[nw/2].set_ylabel(r"$-r^2 \xi(r)\mathrm{[h^{-2}Mpc^2]}$")
else :
    ax[nw/2].set_ylabel(r"$r^2 \xi(r)\mathrm{[h^{-2}Mpc^2]}$")

colors=["b","r","g","k"]


for w,wedge in zip(range(nw),wedges) :
    print "plotting mu",wedge
    
    xidatav=[]
    
    
    first=True
    for d,c,color in zip(data,cov,colors) :
	if args.no_ivar_weight: 
            r,xidata,xierr,wedge_cov=compute_wedge(d,c,murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin)      
	else: 
            r,xidata,xierr,wedge_cov=compute_wedge_with_ivar(d,c,murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin) 
	scale=r**args.rpower
        
        if args.flip :
            ax[w].errorbar(r,-scale*xidata,scale*xierr,fmt="o",color=color)
        else :
            ax[w].errorbar(r,scale*xidata,scale*xierr,fmt="o",color=color)
	ax[w].grid(b=True)
        
        xidatav.append(xidata)
        first=False

    for i in range(len(models)) :
        model=models[i]
        color=colors[i]
        if len(cov)>i :
                c=cov[i]
        else :
                c=cov[0]
	if args.no_ivar_weight: r,ximod,junk,junk=compute_wedge(model,c,murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin)
	else: r,ximod,junk,junk=compute_wedge_with_ivar(model,c,murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin)
        if args.flip :
            ax[w].plot(r,-scale*ximod,"-",color=color,linewidth=2)
        else :
            ax[w].plot(r,scale*ximod,"-",color=color,linewidth=2)
    
    if args.chi2 and len(data)==1 and len(models)==0 :
        weight=np.linalg.inv(wedge_cov)
        res=xidata
        chi2=np.inner(res,weight.dot(res))
        ndata=res.size
        print "(data-0) chi2/ndata=%f/%d=%f"%(chi2,ndata,chi2/ndata)
        
    if args.chi2 and len(data)==1 and len(models)==1 :
        weight=np.linalg.inv(wedge_cov)
        res=xidata-ximod
        chi2=np.inner(res,weight.dot(res))
        ndata=res.size
        print "(data-model) chi2/ndata=%f/%d=%f"%(chi2,ndata,chi2/ndata)
    
    if args.chi2 and len(data)==2 :
        weight=np.linalg.inv(wedge_cov)
        res=xidatav[0]-xidatav[1]
        chi2=np.inner(res,weight.dot(res))
        ndata=res.size
        print "(data0-data1) chi2/ndata=%f/%d=%f"%(chi2,ndata,chi2/ndata)
    
    ax[w].set_title(r"$%2.2f < \mu < %2.2f$"%(wedges[w][0],wedges[w][1]))


plt.xlabel(r"$r\mathrm{[h^{-1}Mpc]}$")

if not args.noshow: plt.show()

if args.out != None:
	f.savefig(args.out+".png",bbox_inches="tight")

