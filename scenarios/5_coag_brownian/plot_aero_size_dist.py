from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

# file name
i_hr = 25
fn = "./out/example_0001_%08i.nc" %i_hr
# load file
f = netcdf.netcdf_file(fn, 'r', mmap=False)

# number conc. of the population
num_conc_per_particle = f.variables["aero_num_conc"][:]

# mass conc. of the population
mass_per_particle = f.variables["aero_particle_mass"][:].sum(axis=0)

# get particle volume and diameter
aero_density = f.variables["aero_density"][:].reshape(20,-1)
aero_particle_mass = f.variables["aero_particle_mass"][:]
aero_volume_per_particle = (aero_particle_mass/aero_density).sum(axis=0)
aero_diameter = np.cbrt(aero_volume_per_particle*6.0/np.pi)

# ## number distribution
# setup the bins range
bins = np.logspace(-8,-6,3*20+1)
plt.hist(aero_diameter, bins=bins, weights=num_conc_per_particle/np.log10(bins[1]/bins[0]))

plt.xscale('log')
plt.xlabel('diameter $D_p$ [m]')
plt.ylabel(r'$dN/d \log D_p$ [# m$^{-3}$]')
plt.savefig('out/size_dist_%03i.pdf' % i_hr)
