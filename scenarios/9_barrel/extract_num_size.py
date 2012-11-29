import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import os
import numpy
import mpl_helper
import matplotlib.pyplot as plt

file = "out_as_0_3LPM_0925/barrel_wc_0001_00000001.nc"
ncf = scipy.io.netcdf.netcdf_file(file, 'r')
particles = partmc.aero_particle_array_t(ncf)
env_state = partmc.env_state_t(ncf)
ncf.close()

dry_diameters = particles.dry_diameters()
x_values = dry_diameters
x_grid = partmc.log_grid(min=1e-8, max=1e-6, n_bin=100)

dist = partmc.histogram_1d(x_values, x_grid, weighted=True, weights=particles.num_concs)

ref_data = numpy.loadtxt("ref_aero_size_num_regrid.txt")

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.semilogx(x_grid.centers(),dist,marker='o',markersize=3.5,color='k')
axes.semilogx(ref_data[:,0],ref_data[:,1],linestyle='-',linewidth=1,color='r')
#axes.plot(x_grid.centers(), dist)
axes.set_title("")
axes.set_xlabel("Dry diameter (m)")
axes.set_ylabel("Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.legend(('PartMC','Barrel'))
filename_out = "aero_num_size.pdf"
figure.savefig(filename_out)
#str_command = "open "+filename_out
#os.system(str_command)
