import partmc
import scipy.io
import numpy
import mpl_helper
import matplotlib.pyplot as plt

file = "barrel_wc_0001_00000001.nc"
ncf = scipy.io.netcdf.netcdf_file(file, 'r')
particles = partmc.aero_particle_array_t(ncf)
env_state = partmc.env_state_t(ncf)
ncf.close()

dry_diameters = particles.dry_diameters()
x_values = dry_diameters
x_grid = partmc.log_grid(min=1e-8, max=1e-6, n_bin=100)

dist = partmc.histogram_1d(x_values, x_grid, weighted=True, weights=particles.num_concs)

(figure, axes) = mpl_helper.make_fig(colorbar=False)
title = ""
xlabel = ""
ylabel = ""
plt.semilogx(x_grid.centers(), dist)
plt.savefig('aero_num_size.pdf')

f1 = 'bin_centers.txt'
numpy.savetxt(f1,x_grid.centers())
f2 = 'num_conc.txt'
numpy.savetxt(f2,dist)

