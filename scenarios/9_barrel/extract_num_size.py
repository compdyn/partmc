import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import os
import numpy
import mpl_helper
import matplotlib.pyplot as plt

file = "out/barrel_wc_0001_00000049.nc"
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
plt.plot(x_grid.centers(), dist)
filename_out = "aero_num_size.pdf"
plt.savefig(filename_out)
str_command = "open "+filename_out
os.system(str_command)
