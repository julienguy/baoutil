#!/usr/bin/env python 
import h5py 
import argparse
from scipy.stats import chi2

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required=True, help = 'h5 input file')
args = parser.parse_args()

f = h5py.File(args.infile)

x=f['best fit'].attrs
print
for k in x :
    print(k)

for k in ['list of fixed pars','list of free pars'] :
    print(k,":")
    for i in x[k]:
        try : 
            if k == 'list of free pars' :
                print("%s = %2.4f +/- %2.4f"%(i, x[i][0], x[i][1]))
            else :
                print("%s = %2.4f"%(i, x[i][0]))
        except KeyError:
            continue

fit_attr = ['fval','ndata','npar']
fval  = x['fval']
ndata = x['ndata']
npar  = x['npar']
df = ndata-npar
p = 1 - chi2.cdf(fval,df)

print
print("chi2 = %2.2f, ndata = %i, npar = %i"%(fval,ndata,npar))
print("chi2/DOF = %2.2f/(%i-%i) = %5.3f"%(fval,ndata,npar,fval/(ndata-npar)))
print("p = %2.3f"%(p))


    

