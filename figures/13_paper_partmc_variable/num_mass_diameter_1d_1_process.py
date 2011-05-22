#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
sys.path.append("../../tool")
import partmc
import config

i_loop_max = config.i_loop_max

def make_plot(in_files, f1, f2, f3, f4, f5):
    x_axis = partmc.log_grid(min=1e-10,max=1e-4,n_bin=100)
    x_centers = x_axis.centers()
    counter = 0
    hist_array_num = np.zeros([len(x_centers),config.i_loop_max])
    hist_array_mass = np.zeros([len(x_centers),config.i_loop_max])
    for file in in_files:
        ncf = scipy.io.netcdf.netcdf_file(config.netcdf_dir+'/'+file, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        ncf.close()

        dry_diameters = particles.dry_diameters()
        hist = partmc.histogram_1d(dry_diameters, x_axis, weights = 1 / particles.comp_vols)
        hist_array_num[:,counter] = hist
        hist = partmc.histogram_1d(dry_diameters, x_axis, weights = particles.masses(exclude=["H2O"])  / particles.comp_vols)
        hist_array_mass[:,counter] = hist

        counter = counter+1

    hist_array_gav_num = np.exp(np.average(np.log(hist_array_num),axis = 1))
    hist_array_gstd_num = np.exp(np.std(np.log(hist_array_num), axis = 1))
    e_bar_top_num = hist_array_gav_num * hist_array_gstd_num
    e_bar_bottom_num = hist_array_gav_num / hist_array_gstd_num
    e_bars_num = np.vstack((hist_array_gav_num - e_bar_bottom_num, e_bar_top_num - hist_array_gav_num))

    hist_array_gav_mass = np.exp(np.average(np.log(hist_array_mass),axis = 1))
    hist_array_gstd_mass = np.exp(np.std(np.log(hist_array_mass), axis = 1))
    e_bar_top_mass = hist_array_gav_mass * hist_array_gstd_mass
    e_bar_bottom_mass = hist_array_gav_mass / hist_array_gstd_mass
    e_bars_mass = np.vstack((hist_array_gav_mass - e_bar_bottom_mass, e_bar_top_mass - hist_array_gav_mass))

    np.savetxt(f1, x_axis.centers())
    np.savetxt(f2, hist_array_gav_num)
    np.savetxt(f3, e_bars_num)
    np.savetxt(f4, hist_array_gav_mass)
    np.savetxt(f5, e_bars_mass)

#used for urban plume scenario data
#for case in ["10K_wei+1", "10K_wei-1", "10K_wei-4" ]:
#    for hour in range(12, 13):
#        files = []
#        for i_loop in range (0, config.i_loop_max):
#            filename_in = "urban_plume_wc_%s_0%03d_000000%02d.nc" % (case, (i_loop+1), hour)
#            files.append(filename_in)
#        f1 = "data/1d_%s_%02d_x_values.txt" % (case, hour)
#        f2 = "data/1d_%s_%02d_hist_array_gav_num.txt" % (case, hour)
#        f3 = "data/1d_%s_%02d_e_bars_num.txt" % (case, hour)
#        f4 = "data/1d_%s_%02d_hist_array_gav_mass.txt" % (case, hour)
#        f5 = "data/1d_%s_%02d_e_bars_mass.txt" % (case, hour)
#        print f1
#        make_plot(files, f1, f2, f3, f4, f5)

for hour in range(12, 13):
    files = []
    for i_loop in range (0, config.i_loop_max):
        filename_in = "urban_plume_wc_%s_0%03d_000000%02d.nc" % (case, (i_loop+1), hour)
        files.append(filename_in)
    f1 = "data/1d_%s_%02d_x_values.txt" % (case, hour)
    f2 = "data/1d_%s_%02d_hist_array_gav_num.txt" % (case, hour)
    f3 = "data/1d_%s_%02d_e_bars_num.txt" % (case, hour)
    f4 = "data/1d_%s_%02d_hist_array_gav_mass.txt" % (case, hour)
    f5 = "data/1d_%s_%02d_e_bars_mass.txt" % (case, hour)
    print f1
    make_plot(files, f1, f2, f3, f4, f5)




    



