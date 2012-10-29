#!/usr/bin/env python

import sys
sys.path.append('../../tool/')
import partmc
import scipy.io
import numpy

prefix = "out/ship_plume_nummass"

grid = partmc.log_grid(1e-9, 1e-3, 300)

n_repeat = 10
n_index = 11

avg_data = numpy.zeros((grid.n_bins, n_index + 1))
for i_repeat in range(1, n_repeat + 1):
    for i_index in range(1, n_index + 1):
        filename = "%s_%04d_%08d.nc" % (prefix, i_repeat, i_index)
        print filename
        ncf = scipy.io.netcdf_file(filename)
        p = partmc.aero_particle_array_t(ncf)
        ncf.close()
        d = p.dry_diameters()
        num_conc = partmc.histogram_1d(d, grid, 1 / p.comp_vols)
        data = numpy.zeros((grid.n_bins, 2))
        data[:,0] = grid.centers()
        data[:,1] = num_conc
        out_filename = '%s_%04d_%08d_num.txt' % (prefix, i_repeat, i_index)
        numpy.savetxt(out_filename, data)
        avg_data[:, i_index + 1] += num_conc

avg_data /= n_repeat
avg_data[:,0] = grid.centers()
out_filename = '%s_avg_num.txt' % (prefix, i_repeat, i_index)
numpy.savetxt(out_filename, avg_data)
