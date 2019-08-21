#!/usr/bin/env python 
import h5py 
import argparse
from scipy.stats import chi2

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required=True, help = 'h5 input file')
args = parser.parse_args()

f = h5py.File(args.infile,mode="r")

x=f['best fit'].attrs
#print
#for k in x :
#    print(k)

for k in ['list of fixed pars','list of free pars'] :
    print("\n"+k)
    print("-------------------------------------------------")
    for i in x[k]:
        try : 
            if k == 'list of free pars' :
                print("%s = %2.4f +/- %2.4f"%(i, x[i][0], x[i][1]))
                param=i.decode("utf-8")
                if param.find("bias_eta")>=0 :
                    what=param.replace("bias_eta_","")
                    beta=x["beta_{}".format(what)][0]
                    growthrate=x["growth_rate"][0]
                    bias=x[i][0]/beta*growthrate
                    bias_err=x[i][1]/beta*growthrate
                    print("bias_%s = %2.4f +/- %2.4f (using f=%2.4f beta=%2.4f)"%(what, bias, bias_err,growthrate, beta))
            else :
                print("%s = %2.4f"%(i, x[i][0]))
                param=i.decode("utf-8")
                if param.find("bias_eta")>=0 :
                    what=param.replace("bias_eta_","")
                    beta=x["beta_{}".format(what)][0]
                    betaerr=x["beta_{}".format(what)][1]
                    growthrate=x["growth_rate"][0]
                    bias=x[i][0]/beta*growthrate
                    bias_err=x[i][0]*betaerr/beta**2*growthrate
                    print("bias_%s = %2.4f +/- %2.4f (using f=%2.4f beta=%2.4f)"%(what, bias, bias_err,growthrate, beta))
        except KeyError:
            continue

fit_attr = ['fval','ndata','npar']
fval  = x['fval']
ndata = x['ndata']
npar  = x['npar']
df = ndata-npar
p = 1 - chi2.cdf(fval,df)

print("\nchi2")
print("-------------------------------------------------")    
print("chi2 = %2.2f, ndata = %i, npar = %i"%(fval,ndata,npar))
print("chi2/DOF = %2.2f/(%i-%i) = %5.3f"%(fval,ndata,npar,fval/(ndata-npar)))
print("p = %2.3f"%(p))


    

