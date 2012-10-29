#!/usr/bin/env python

import sys
sys.path.append('../../../tool/')
import partmc
import scipy.io
import numpy

grid = partmc.log_grid(1e-9, 1e-3, 300)

for (t, filename, key) in partmc.get_time_filename_list('.', r'ship_plume_wc_0001_(.*).nc'):
    print filename
    ncf = scipy.io.netcdf_file(filename)
    p = partmc.aero_particle_array_t(ncf)
    ncf.close()
    d = p.diameters()
    w = p.weight_groups
    nc1 = partmc.histogram_1d(d[w == 1], grid) * grid.grid_size(0)
    nc2 = partmc.histogram_1d(d[w == 2], grid) * grid.grid_size(0)
    c1 = p.comp_vols[w == 1]
    c2 = p.comp_vols[w == 2]
    n1 = partmc.histogram_1d(d[w == 1], grid, 1 / c1)
    n2 = partmc.histogram_1d(d[w == 2], grid, 1 / c2)
    data = numpy.zeros((len(n1), 5))
    data[:,0] = grid.centers()
    data[:,1] = nc1
    data[:,2] = nc2
    data[:,3] = n1
    data[:,4] = n2
    out_filename = 'comp_num_%s.txt' % key
    numpy.savetxt(out_filename, data)
