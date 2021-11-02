from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt

# ## create a function and plot the time series
def get_num_conc(fn):
    with netcdf.netcdf_file(fn, 'r') as f:
        return f.variables["aero_num_conc"].data.sum()

from glob import glob
fn_ls = sorted(glob('./out/example_*.nc'))

num_conc_ls = []
for fn in fn_ls:
    num_conc_ls.append(get_num_conc(fn))
    
plt.plot(num_conc_ls)
plt.xlabel("hour")
plt.ylabel(r"number conc [# m$^{-3}$]")
plt.savefig('out/num_conc.pdf')
