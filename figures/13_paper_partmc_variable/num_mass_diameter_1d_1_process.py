#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
sys.path.append("../../tool")
import partmc
import config

i_loop_max = config.i_loop_max

def make_plot(in_files, f1, f2, f3):
    x_axis = partmc.log_grid(min=1e-10,max=1e-4,n_bin=100)
    x_centers = x_axis.centers()
    counter = 0
    hist_array = np.zeros([len(x_centers),config.i_loop_max])
    for file in in_files:
        ncf = scipy.io.netcdf.netcdf_file(config.netcdf_dir+file, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        ncf.close()

        dry_diameters = particles.dry_diameters()
        hist = partmc.histogram_1d(dry_diameters, x_axis, weights = 1 / particles.comp_vols)
        hist_array[:,counter] = hist
        counter = counter+1

    hist_array_gav = np.exp(np.average(np.log(hist_array),axis = 1))
    hist_array_gstd = np.exp(np.std(np.log(hist_array), axis = 1))
    e_bar_top = hist_array_gav * hist_array_gstd
    e_bar_bottom = hist_array_gav / hist_array_gstd
    e_bars = np.vstack((hist_array_gav - e_bar_bottom, e_bar_top - hist_array_gav))

    np.savetxt(f1, x_axis.centers())
    np.savetxt(f2, hist_array_gav)
    np.savetxt(f3, e_bars)

for case in ["10K_wei+1", "10K_wei-1", "10K_wei-4" ]:
    for hour in range(12, 13):
        files = []
        for i_loop in range (0, config.i_loop_max):
            filename_in = "urban_plume_wc_%s_0%03d_000000%02d.nc" % (case, (i_loop+1), hour)
            files.append(filename_in)
        f1 = "data/1d_%s_%02d_x_values.txt" % (case, hour)
        f2 = "data/1d_%s_%02d_hist_array_gav.txt" % (case, hour)
        f3 = "data/1d_%s_%02d_e_bars.txt" % (case, hour)
        print f1
        make_plot(files, f1, f2, f3)




    



