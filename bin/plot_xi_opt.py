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
from baoutil.wedge import optprofile

plt.rcParams["font.family"]="serif"
plt.rcParams["font.size"]=16.0


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-d','--data', type = str, default = None, required=True, nargs='*',
                        help = 'baofit data')
parser.add_argument('--rrange', type = str, default = "10:180", required=False,
                        help = 'r range for wedge of the form "rmin:rmax')
parser.add_argument('--rbin', type = float, default = 4.0, required=False,
                        help = 'r bin size')
parser.add_argument('--rpmin', type = float, default = 0, required=False,
                        help = 'min r_parallel')
parser.add_argument('--res', type = str, default = None, required=False, nargs='*',
                    help = 'baofit residuals file to plot model')

parser.add_argument('--model', type = str, default = None, required=False, nargs='*',
                    help = 'model file')

parser.add_argument('--out', type = str, default = None, required=False,
		                        help = 'output prefix')
parser.add_argument('--chi2', action="store_true",
		help = 'compute chi2 of wedges (if only one data set and model)')
parser.add_argument('--out-txt', type = str, default = None, required=False,
		                        help = 'output text file')

parser.add_argument('--noshow', action="store_true",
		            help = 'prevent the figure window from displaying')
parser.add_argument('--rpower', type = int, default = 2, required=False,
                        help = 'r power for display')
parser.add_argument('--flip', action="store_true",
		            help = 'flip plot (useful for Lya-QSO cross-corr)')
parser.add_argument('--beta', type=float, default=1.5, 
		            help = 'beta value for Kaiser weight in profile')

#parser.add_argument('--ivar_weight', action="store_true",
# help = 'use inverse variance to combine the bins')
parser.add_argument('--no_ivar_weight', action="store_true",
                    help = 'do not use inverse variance to combine the bins')


args = parser.parse_args()


data=[]
cov=[]
for d1 in args.data :
        if d1.find(".fits")>=0 :
                trp,trt,d,c=read_baofit_fits(d1)
                data.append(d)
                cov.append(c)
        else :
            d=read_baofit_data(d1)
            data.append(d)
            if d1.find(".data")>=0 :
                cov_filename=d1.replace(".data",".cov")
                if not os.path.isfile(cov_filename) :
                        cov_filename=string.replace(d1,".data","-cov.fits")
                if os.path.isfile(cov_filename) :
                        c  = read_baofit_cov(cov_filename,n2d=d.size,convert=True)
                        cov.append(c)
                else :
                    print("warning, cannot find covariance of ",d1)
                    cov.append("None")    
            else :
                print("warning, cannot guess covariance of ",d1)
                cov.append("None")


models=[]
if args.res is not None :
    for i,res in enumerate(args.res) :
        if i<len(data) :
            d=data[i]
        else :
            d=data[0]
        mod = read_baofit_model(res,n2d=d.size)
        if d.size != mod.size :
            print("error data and model don't have same size")
            sys.exit(12)
        models.append(mod)
    




if args.rrange :
    try :
        vals=args.rrange.split(":")
        if len(vals)!=2 :
            print("incorrect format for r range '%s', expect rmin:rmax"%args.rrange)
            sys.exit(12)
        rmin=float(vals[0])
        rmax=float(vals[1])
        rrange=[rmin,rmax]
    except ValueError as e:
        print(e)
        print("incorrect format for r range '%s', expect rmin:rmax"%args.rrange)
        sys.exit(12)
else :
    rrange=[10,180]

f = plt.figure()
ax = plt.subplot(1,1,1)
ax.set_ylabel(r"$r^2 \xi(r)\mathrm{[h^{-2}Mpc^2]}$")

data_colors=["k","b","r","g","k"]
model_colors=["b","b","b","gray"]
model_alphas=[1,0.6,0.4,0.4]
xidatav=[]

rout=None
yout=None
eout=None
mout=None

    
first=True
for d,c,color in zip(data,cov,data_colors) :
        r,xidata,xierr,wedge_cov=optprofile(d,c,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin,beta=args.beta) 
        scale=r**args.rpower
        ax.errorbar(r,scale*xidata,scale*xierr,fmt="o",color=color)
        ax.grid(b=True)
        if args.out_txt :
                rout=r
                yout=xidata
                eout=xierr
                

for i in range(len(models)) :
        model=models[i]
        color=model_colors[i]
        alpha=model_alphas[i]
        if len(cov)>i :
                c=cov[i]
        else :
                c=cov[0]
        r,ximod,junk,junk=optprofile(model,c,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin,beta=args.beta)
        ax.plot(r[ximod!=0],(scale*ximod)[ximod!=0],"-",color=color,linewidth=2,alpha=alpha)
        if args.out_txt :
                mout = ximod

ax.set_xlabel(r"$r\mathrm{[h^{-1}Mpc]}$")

if not args.noshow: plt.show()

if args.out != None:
	f.savefig(args.out+".png",bbox_inches="tight")

if args.out_txt is not None :
        if mout is None :
                tmp=np.array([rout,yout,eout])
        else :
                tmp=np.array([rout,yout,eout,mout])
        np.savetxt(args.out_txt,tmp.T)
        print("wrote",args.out_txt)
