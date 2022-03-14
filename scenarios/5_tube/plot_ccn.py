import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import matplotlib as mpl
mpl.rcParams['font.size'] = 12

#Read files
out_filename = 'out_as2suc/urban_plume_ccn_20_lin.pdf'
filename = 'out_as2suc/urban_plume_00000020_process.nc'
ncf  = scipy.io.netcdf_file(filename, mmap=False)
sc_2d_dist   = ncf.variables["diam_sc_dist"].data * 1e-6
wet_num_dist = ncf.variables["num_dist"].data * 1e-6
sc_grid      = ncf.variables["sc"].data * 100 # Supersaturation level in %
diam         = ncf.variables["diam"].data*1e6 # diameter in um

print("sc_grid")
print(sc_grid)

# Caculate the cumulative sum of each ss level
ccn = np.nancumsum(sc_2d_dist, axis=0)*np.log(sc_grid[1]/sc_grid[0]) 
plt.figure(figsize=(9,8))
plt.subplot(2,1,1)
plt.xscale('log');
plt.plot(diam, wet_num_dist, label='Number distribution',color='k')
for i in (5,10,15,20):
    ss = sc_grid[i]
    plt.xlim(1e-2,1e0)
    plt.plot(diam, ccn[i], label="CCN"+'ss=%.2f%%'%ss + ' distribution',ls='--')
plt.legend(loc=1)
plt.grid()
plt.ylabel(r'Number concentration, $\rm cm^{-3}$')

plt.subplot(2,1,2)
plt.xscale('linear');
for i in (5,10,15,20):
    ss = sc_grid[i]
    plt.xlim(1e-2,1e0)
    plt.plot(diam, ccn[i]/wet_num_dist, label="CCN"+'ss=%.2f%%'%ss)
    plt.legend(loc=4)
plt.xlabel(r'Diameter, $\rm \mu m$')
plt.ylabel(r'$N_{\rm CCN}/N_{\rm CN}$')
plt.grid()
plt.savefig(out_filename)
print(out_filename)

