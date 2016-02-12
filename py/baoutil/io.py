import numpy as np
import astropy.io.fits as fits
import string
import os.path
import sys
from math import sqrt

def read_baofit_data(filename) :
    if filename.find(".data")<0 :
        data_filename="%s.data"%filename
    else :
        data_filename=filename
    print "reading data in %s"%data_filename

    if not os.path.isfile(data_filename) :
        print "error %s doesn't exist"%data_filename
        sys.exit(12)

    vals=np.loadtxt(data_filename).T
    return vals[1]

def read_baofit_model(res_filename,n2d=2500) :
    if not os.path.isfile(res_filename) :
        print "error %s doesn't exist"%res_filename
        sys.exit(12)

    vals=np.loadtxt(res_filename).T
    i=vals[0].astype(int)
    ni=np.max(i)+1
    if ni>n2d :
        n2d=ni
    m=vals[7]
    res=np.zeros((n2d))
    res[i]=m
    return res


def read_baofit_cov(filename,convert=True) :
    print "reading cov in %s"%filename
    
    if filename.find(".fits")>0 :
        return fits.open(filename)[0].data
    
    if filename.find(".cov")>0 :
        cov_filename=filename
        base_filename=string.replace(filename,".cov","")
        fits_filename="%s-cov.fits"%base_filename
    else :
        base_filename=filename
        cov_filename="%s.cov"%base_filename
        fits_filename="%s-cov.fits"%base_filename
    
    if os.path.isfile(fits_filename) :
        print "using %s"%fits_filename
        print "NEED TO CHECK DATE"
        return fits.open(fits_filename)[0].data
    
    if not os.path.isfile(cov_filename) :
        print "error %s doesn't exist"%cov_filename
        sys.exit(12)

    print "reading %s"%cov_filename
    vals=np.loadtxt(cov_filename).T
    
    ii=vals[0].astype(int)
    jj=vals[1].astype(int)
    cc=vals[2]
    
    n=max(np.max(ii),np.max(jj))+1
    cov=np.zeros((n,n))
    
    for i,j,c in zip(ii,jj,cc) :
        cov[i,j]=c
        cov[j,i]=c
    
    if convert :
        fits.writeto(fits_filename,cov)
        print "wrote %s"%fits_filename
    
    return cov

