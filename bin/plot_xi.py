#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

import os.path

from baoutil.io import read_baofit_data,read_baofit_cov,read_baofit_fits,read_baofit_model
from baoutil.wedge import compute_wedge,compute_wedge_with_ivar,block

def getlabel(wedge,sign) :
        if rp_sign>0 :
                return "$%3.2f<\mu<%3.2f$"%(wedge[0],wedge[1])
        else :
                return "$%3.2f<\mu<%3.2f$"%(-wedge[1],-wedge[0])

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

parser.add_argument('--abs', action="store_true",
                    help = 'average as a function of |rp|')


args = parser.parse_args()

rp=None
rt=None
data=[]
cov=[]
for d1 in args.data :
        if d1.find(".fits")>=0 :
                trp,trt,d,c=read_baofit_fits(d1)
                if args.abs :
                        trp=np.abs(trp)
                data.append(d)
                cov.append(c)
                if rp is None :
                        rp=trp
                        rt=trt                        
                else :
                        if np.max(np.abs((rp-trp)))>01 :
                                print("rp values don't match")
                                sys.exit(12)
                        if np.max(np.abs((rt-trt)))>01 :
                                print("rt values don't match")
                                sys.exit(12)
        else :
                d=read_baofit_data(d1)

                if rp is None :
                      print("warning, making up rp and rt values")
                      n2d=d.size
                      n1d=np.sqrt(n2d).astype(int)
                      rstep=4.
                      rt=((np.arange(n2d)%n1d+0.5)*rstep).astype(float)
                      rp=((np.arange(n2d)/n1d+0.5)*rstep).astype(float)
                
                data.append(d)
                if d1.find(".data")>=0 :
                        cov_filename=d1.replace(".data",".cov")
                        if not os.path.isfile(cov_filename) :
                                cov_filename=d1.replace(".data","-cov.fits")
                        if os.path.isfile(cov_filename) :
                                c  = read_baofit_cov(cov_filename,n2d=d.size,convert=True)
                                cov.append(c)
                        else :
                                print("warning, cannot find covariance of ",d1)
                                cov.append(np.eye(d.size)*1e-12)    
                else :
                        print("warning, cannot guess covariance of ",d1)
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
            print("error data and model don't have same size")
            sys.exit(12)
        models.append(mod)
    

print("rp range = %f %f"%(np.min(rp),np.max(rp)))
print("rt range = %f %f"%(np.min(rt),np.max(rt)))




if args.mu :
    wedges=[]
    try :
        wedge_strings=args.mu.split(",")
        for wedge_string in wedge_strings :
            print(wedge_string)
            vals=wedge_string.split(":")
            if len(vals)!=2 :
                print("incorrect format for mu range '%s', expect mumin:mumax"%args.mu)
                sys.exit(12)
            mumin=float(vals[0])
            mumax=float(vals[1])
            wedges.append([mumin,mumax])
    except ValueError as e:
        print(e)
        print("incorrect format for mu range '%s', expect mumin:mumax"%args.mu)
        sys.exit(12)
else :
    wedges= [[0.8,1.0],[0.5,0.8],[0.0,0.5]]
    wedges= [[0.95,1.0],[0.8,0.95],[0.5,0.8],[0.0,0.5]]

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


nw=len(wedges)

plt.figure()
ncols=int(np.sqrt(nw))
nrows=nw//ncols
if ncols*nrows < nw : nrows += 1
ax=[]
for index in range(1,nw+1) :
        ax.append(plt.subplot(nrows,ncols,index))


if args.flip :
    ax[nw//2].set_ylabel(r"$-r^2 \xi(r)\mathrm{[h^{-2}Mpc^2]}$")
else :
    ax[0].set_ylabel(r"$r^2 \xi(r)\mathrm{[h^{-2}Mpc^2]}$")
    ax[nw//2].set_ylabel(r"$r^2 \xi(r)\mathrm{[h^{-2}Mpc^2]}$")

colors=["b","r","g","k","gray","purple"]

abs_rp=True
has_neg_rp= np.sum(rp<0)>0
rp_signs=[]
if has_neg_rp and abs_rp :        
        rp_signs=[1,-1]
else :
        rp_signs=[1]

xidata_array={}
ximod_array={}
cov_array={}
color_array={}

for w,wedge in zip(range(nw),wedges) :
    print("plotting mu",wedge)
        
    first=True
    color_index=0
    for d,c in zip(data,cov) :

        for rp_sign in rp_signs :
                label = getlabel(wedge,rp_sign)
                subsample=np.where(rp_sign*rp>=0)[0]
                        
                color=colors[color_index]
                color_index+=1
                
                if args.no_ivar_weight : 
                        r,xidata,xierr,wedge_cov=compute_wedge(rp[subsample],rt[subsample],d[subsample],block(c,subsample),murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin)      
                else : 
                        r,xidata,xierr,wedge_cov=compute_wedge_with_ivar(rp[subsample],rt[subsample],d[subsample],block(c,subsample),murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin)

                
                xidata_array[label] = xidata
                cov_array[label] = wedge_cov
                color_array[label] = color
                
                scale=r**args.rpower

                if args.flip :
                    ax[w].errorbar(r,-scale*xidata,scale*xierr,fmt="o",color=color,label=label)
                else :
                    ax[w].errorbar(r,scale*xidata,scale*xierr,fmt="o",color=color,label=label)
                ax[w].grid(b=True)
                ax[w].legend(fontsize="small",numpoints=1,loc="lower left")
                first=False
        
        for i in range(len(models)) :      
                model=models[i]
                if len(cov)>i :
                        c=cov[i]
                else :
                        c=cov[0]
                for rp_sign in rp_signs :
                        label = getlabel(wedge,rp_sign)
                        subsample=np.where(rp_sign*rp>=0)[0]
                        color=color_array[label]
                        
                        if args.no_ivar_weight: 
                                r,ximod,junk,junk=compute_wedge(rp[subsample],rt[subsample],model[subsample],block(c,subsample),murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin)
                        else: 
                                r,ximod,junk,junk=compute_wedge_with_ivar(rp[subsample],rt[subsample],model[subsample],block(c,subsample),murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin)
                        if args.flip :
                                ax[w].plot(r,-scale*ximod,"-",color=color,linewidth=2)
                        else :
                                ax[w].plot(r,scale*ximod,"-",color=color,linewidth=2)

                        ximod_array[label] = ximod
                        
    if args.chi2 and len(data)==1 and len(models)==0 :
           
            weight=np.linalg.inv(wedge_cov)
            res=xidata
            chi2=np.inner(res,weight.dot(res))
            ndata=res.size
            print("(data-0) chi2/ndata=%f/%d=%f"%(chi2,ndata,chi2/ndata))
        
    if args.chi2 and len(data)==1 and len(models)==1 :
            
            weight=np.linalg.inv(wedge_cov)
            res=xidata-ximod
            chi2=np.inner(res,weight.dot(res))
            ndata=res.size
            print("(data-model) chi2/ndata=%f/%d=%f"%(chi2,ndata,chi2/ndata))
    
    if args.chi2 and len(data)==2 :
                    
            weight=np.linalg.inv(wedge_cov)
            res=xidatav[0]-xidatav[1]
            chi2=np.inner(res,weight.dot(res))
            ndata=res.size
            print("(data0-data1) chi2/ndata=%f/%d=%f"%(chi2,ndata,chi2/ndata))

    #ax[w].set_title(r"$%2.2f < \mu < %2.2f$"%(wedges[w][0],wedges[w][1]))




for a in ax[nw-2:nw] :
        a.set_xlabel(r"$r\mathrm{[h^{-1}Mpc]}$")   
if not args.noshow: plt.show()

if args.out != None:
	f.savefig(args.out+".png",bbox_inches="tight")

