#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import string

from baoutil.io import read_baofit_data
from baoutil.io import read_baofit_cov
from baoutil.io import read_baofit_model
from baoutil.wedge import compute_wedge2

plt.rcParams["font.family"]="serif"
plt.rcParams["font.size"]=16.0


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-d','--data', type = str, default = None, required=True, nargs='*', action='append',
                        help = 'baofit data')
parser.add_argument('-c','--cov', type = str, default = None, required=False,
                        help = 'baofit cov (default is guessed from data)')
parser.add_argument('--mu', type = str, default = None, required=False,
                        help = 'mu range for wedge of the form "mumin:mumax,mumin:mumax ...')
parser.add_argument('--np', type = float, default = None, required=False,
                        help = 'n_p, default : 50')
parser.add_argument('--nt', type = float, default = None, required=False,
                        help = 'n_t, default : 50')
parser.add_argument('--rpmin', type = float, default = 10, required=False,
                        help = 'rpmin')
parser.add_argument('--rpmax', type = float, default = 180, required=False,
                        help = 'rpmax')
parser.add_argument('--rbin', type = float, default = 4.0, required=False,
                        help = 'r bin size')
parser.add_argument('--res', type = str, default = None, required=False, nargs='*', action='append', 
                    help = 'baofit residuals file to plot model')
parser.add_argument('--out', type = str, default = None, required=False,
		                        help = 'output prefix')
parser.add_argument('--chi2', action="store_true",
		help = 'compute chi2 of wedges (if only one data set and model)')
parser.add_argument('--noshow', action="store_true",
		            help = 'prevent the figure window from displaying')
parser.add_argument('--flip', action="store_true",
		            help = 'flip plot (useful for Lya-QSO cross-corr)')


args = parser.parse_args()

cov_filename=args.cov



data=[]
name=[]
cov = []
for d2 in args.data : # If more than one file 
    for d1 in d2 :
       data.append(read_baofit_data(d1))
       name.append(d2)

print "number of file =", len(name)

for i in range(len(name)):
    cov_filename=string.replace(args.data[i][0],".data",".cov")
    cov.append(read_baofit_cov(cov_filename,n2d=data[i].size,convert=True))
    print "for", name[i][0], ", ",len(data[i]), "points"

if args.np is not None :
    n_p = args.np 
else : 
    n_p = 50

if args.nt is not None :
    n_t = args.nt 
else : 
    n_t = 50


if args.rpmin :
    rpmin = args.rpmin
else :
    rpmin = 10.

print "rpmin =", rpmin 


if args.rpmax :
    rpmax = args.rpmax
else :
    rpmax = 180.

print "rpmax =", rpmax 



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

data_colors=["b","r","g","k"]
model_colors=["r","k","k","k"]

for w,wedge in zip(range(nw),wedges) :
    print "plotting mu",wedge
    
    xidatav=[]
    
            
    for lab,d,c,i in zip(name,data,data_colors, range(len(name))) :
        rp,xidata,xierr,wedge_cov=compute_wedge2(d,n_t,n_p,cov[i],murange=wedge,rpar_min=rpmin, rpar_max=rpmax,rbin=args.rbin)      
        if len(data[i]) == 5000: 
            rp = rp-200
        scale=1
        if args.flip :
            ax[w].errorbar(rp,-scale*xidata,scale*xierr,fmt="o-",color=c, label=lab[0])
        else :
            ax[w].errorbar(rp,scale*xidata,scale*xierr,fmt="o-",color=c, label=lab[0])
	ax[w].grid(b=True)
        
        xidatav.append(xidata)
        first=False
    
    for lab,model,c, i in zip(name, models, model_colors, range(len(name)))  :
	r,ximod,junk,junk=compute_wedge2(d,n_t,n_p,cov[i],murange=wedge,rpar_min=rpmin, rpar_max=rpmax,rbin=args.rbin)
        if args.flip :
            ax[w].plot(r,-scale*ximod,"-",color=c, label=lab[0])
        else :
            ax[w].plot(r,scale*ximod,"-",color=c, label=lab[0])
    
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


plt.xlabel(r"$r_{||}\mathrm{[h^{-1}Mpc]}$")

if not args.noshow: 
    plt.legend()
    plt.show() 

if args.out != None:
	f.savefig(args.out+".png",bbox_inches="tight")

