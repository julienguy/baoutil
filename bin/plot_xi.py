#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

import os.path

from baoutil.io import read_baofit_data,read_baofit_cov,read_baofit_fits,read_baofit_model
from baoutil.wedge import compute_wedge,block

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
parser.add_argument('--rpmin', type = float, default = None, required=False,
                        help = 'min r_parallel')
parser.add_argument('--rpmax', type = float, default = None, required=False,
                        help = 'max r_parallel')

parser.add_argument('--res', type = str, default = None, required=False, nargs='*',
                    help = 'baofit residuals file to plot model')

parser.add_argument('--model', type = str, default = None, required=False, nargs='*',
                    help = 'model file')

parser.add_argument('--out-figure', type = str, default = None, required=False,
		                        help = 'output prefix')

parser.add_argument('--out-txt', type = str, default = None, required=False,
		                        help = 'output text file')
parser.add_argument('--chi2', action="store_true",
		help = 'compute chi2 of wedges (if only one data set and model)')

parser.add_argument('--noshow', action="store_true",
		            help = 'prevent the figure window from displaying')
parser.add_argument('--rpower', type = int, default = 2, required=False,
                        help = 'r power for display')
parser.add_argument('--flip', action="store_true",
		            help = 'flip plot (useful for Lya-QSO cross-corr)')
parser.add_argument('--single-plot', action="store_true",
		            help = 'all wedges on same plot')
parser.add_argument('--legend-loc', type=str, default="lower right" , required = False,
		            help = 'legend location')
parser.add_argument('--title', type=str, default=None , required = False,
		            help = 'title of first subplot')
parser.add_argument('--labels', type=str, default=None , required = False, nargs="*",
		            help = 'labels')

#parser.add_argument('--ivar_weight', action="store_true",
# help = 'use inverse variance to combine the bins')
parser.add_argument('--no_ivar_weight', action="store_true",
                    help = 'do not use inverse variance to combine the bins')

parser.add_argument('--abs', action="store_true",
                    help = 'average as a function of |rp|')
parser.add_argument('--what', type=str, default="LYA(LYA)xLYA(LYA)",
                    required=False,
                    help = 'key to find model in h5 file')
parser.add_argument('--beta', type=float, default=0, 
		            help = 'beta value for Kaiser weight in profile')
parser.add_argument('--same-color', action='store_true', 
		            help = 'force same color for model and data')

args = parser.parse_args()

rp=None
rt=None
data=[]
cov=[]
for dataindex,d1 in enumerate(args.data) :
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
                        continue
                        if np.max(np.abs((rp-trp)))>4. :
                                print("rp values don't match")
                                sys.exit(12)
                        if np.max(np.abs((rt-trt)))>4. :
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
        mod = read_baofit_model(res,n2d=d.size,what=args.what)
        if d.size != mod.size :
            print("error data and model don't have same size")
            sys.exit(12)
        models.append(mod)

if args.model is not None :
    for i,filename in enumerate(args.model) :
        if i<len(data) :
            d=data[i]
        else :
            d=data[0]
        mod = np.loadtxt(filename)
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
    #wedges= [[0.95,1.0],[0.8,0.95],[0.5,0.8],[0.0,0.5]]

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
ax={}
if args.single_plot :
        ncols=1
        nrows=1
        a = plt.subplot(1,1,1)
        for i in range(nw) :
                ax[i]=a
else :
        ncols=int(np.sqrt(nw))
        nrows=nw//ncols
        if ncols*nrows < nw : nrows += 1
        for index in range(1,nw+1) :
                ax[index-1] = plt.subplot(nrows,ncols,index)


if args.flip :
        sign="-"
else :
        sign=""
if args.rpower==0 : 
        label="$%s\\xi(r)$"%sign
else :
        p=args.rpower
        label="%s$r^%d \\xi(r)\\mathrm{[h^{-%d}Mpc^%d]}$"%(sign,p,p,p)


for i in range(0,nw,ncols) :
        ax[i].set_ylabel(r"%s"%label)


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

color_index=0
other_color_index=0

rout=None
yout=None
eout=None
mout=None

for w,wedge in zip(range(nw),wedges) :
    print("plotting mu",wedge)
        
    first=True
    if not args.single_plot : color_index=0 
    
    for dataindex,d,c in zip(np.arange(len(data)),data,cov) :

        for rp_sign in rp_signs :
                label = getlabel(wedge,rp_sign)
                if args.labels is not None and len(args.labels) == len(data) :
                        full_label = args.labels[dataindex]+" "+label
                else :
                        full_label = label
                subsample=np.where(rp_sign*rp>=0)[0]
                        
                color=colors[color_index]
                color_index+=1
                color_index = color_index%len(colors)
                r,xidata,xierr,wedge_cov=compute_wedge(rp[subsample],rt[subsample],d[subsample],block(c,subsample),murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin,beta=args.beta,rpmax=args.rpmax)

                        
                if args.out_txt :
                        rout=r
                        yout=xidata
                        eout=xierr
                
                        
                
                xidata_array[label] = xidata
                cov_array[label] = wedge_cov
                color_array[label] = color
                
                scale=r**args.rpower

                if args.flip :
                    ax[w].errorbar(r,-scale*xidata,scale*xierr,fmt="o",color=color,label=full_label)
                else :
                    ax[w].errorbar(r,scale*xidata,scale*xierr,fmt="o",color=color,label=full_label)
                ax[w].grid(b=True)
                ax[w].legend(fontsize="small",numpoints=1,loc=args.legend_loc)
                first=False
        if not args.single_plot : other_color_index=0 
        for i in range(len(models)) :      
                model=models[i]
                if len(cov)>i :
                        c=cov[i]
                else :
                        c=cov[0]
                for rp_sign in rp_signs :
                        label = getlabel(wedge,rp_sign)
                        
                        subsample=np.where(rp_sign*rp>=0)[0]
                        if len(data)>1 :
                                color=color_array[label]
                        else :
                                color=colors[color_index]
                                color_index+=1
                                color_index = color_index%len(colors)
                        if args.same_color :
                               color=colors[other_color_index]
                               other_color_index+=1
                               other_color_index = other_color_index%len(colors)  
                        r,ximod,junk,junk=compute_wedge(rp[subsample],rt[subsample],model[subsample],block(c,subsample),murange=wedge,rrange=rrange,rbin=args.rbin,rpmin=args.rpmin,beta=args.beta)
                        if args.flip :
                                ax[w].plot(r,-scale*ximod,"-",color=color,linewidth=2)
                        else :
                                ax[w].plot(r,scale*ximod,"-",color=color,linewidth=2)

                        ximod_array[label] = ximod

                if args.out_txt :
                        mout = ximod

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
            if 1 :
                print("writing this in chi2.fits")
                import astropy.io.fits as pyfits
                h=pyfits.HDUList([pyfits.PrimaryHDU(weight),pyfits.ImageHDU(xidata,name="DATA"),pyfits.ImageHDU(ximod,name="MODEL")])
                h[0].header["EXTNAME"]="WEIGHT"
                h.writeto("chi2.fits",overwrite=True)
                
    if args.chi2 and len(data)==2 :
                    
            weight=np.linalg.inv(wedge_cov)
            res=xidatav[0]-xidatav[1]
            chi2=np.inner(res,weight.dot(res))
            ndata=res.size
            print("(data0-data1) chi2/ndata=%f/%d=%f"%(chi2,ndata,chi2/ndata))

    #ax[w].set_title(r"$%2.2f < \mu < %2.2f$"%(wedges[w][0],wedges[w][1]))


if args.title is not None :
        ax[0].set_title(args.title)

for i in range(nw-ncols,nw) :
        ax[i].set_xlabel(r"$r\mathrm{[h^{-1}Mpc]}$")

if not args.noshow: plt.show()

if args.out_figure != None:
	f.savefig(args.out_figure+".png",bbox_inches="tight")


if args.out_txt is not None :
        if mout is None :
                tmp=np.array([rout,yout,eout])
        else :
                tmp=np.array([rout,yout,eout,mout])
        np.savetxt(args.out_txt,tmp.T)
        print("wrote",args.out_txt)
                
