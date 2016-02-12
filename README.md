# baoutil
python scripts for Lya BAO analysis

usage

```

export PATH=$SOFTDIR/baoutil/bin:${PATH}
export PYTHONPATH=$SOFTDIR/baoutil/py:${PYTHONPATH}

plot_xi.py --data XXX --res XXX_residuals.dat

```

where XXX.data is a baofit input ASCII file of a 2D correlation function 

The script finds the covariance XXX.cov or XXX-cov.fits, and converts the ASCII file to fits if necessary to read it faster next time.


