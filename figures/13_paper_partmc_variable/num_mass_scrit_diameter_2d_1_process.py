#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
sys.path.append("../../tool")
import partmc
import config

i_loop_max = config.i_loop_max

def make_plot(in_files, f1, f2, f3, f4, f5, f6):
    x_axis = partmc.log_grid(min=1e-9,max=1e-5,n_bin=70)
    y_axis = partmc.log_grid(min=1e-3,max=1e2,n_bin=50)
    x_centers = x_axis.centers()
    y_centers = y_axis.centers()
    counter = 0
    hist_array_num = np.zeros([len(x_centers),len(y_centers),config.i_loop_max])
    hist_average_num = np.zeros([len(x_centers), len(y_centers)])
    hist_std_num = np.zeros([len(x_centers), len(y_centers)])
    hist_std_norm_num = np.zeros([len(x_centers), len(y_centers)])

    hist_array_mass = np.zeros([len(x_centers),len(y_centers),config.i_loop_max])
    hist_average_num = np.zeros([len(x_centers), len(y_centers)])
    hist_std_num = np.zeros([len(x_centers), len(y_centers)])
    hist_std_norm_num = np.zeros([len(x_centers), len(y_centers)])
    for file in in_files:
        ncf = scipy.io.netcdf.netcdf_file(config.netcdf_dir+'/'+file, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        s_crit = (particles.critical_rel_humids(env_state) - 1)*100
        dry_mass = particles.masses(exclude = ["H2O"])

        dry_diameters = particles.dry_diameters()
        hist2d = partmc.histogram_2d(dry_diameters, s_crit, x_axis, y_axis, weights = 1 / particles.comp_vols)
        hist_array_num[:,:,counter] = hist2d
        hist2d = partmc.histogram_2d(dry_diameters, s_crit, x_axis, y_axis, weights = particles.masses(include=["BC"])  / particles.comp_vols)
        hist_array_mass[:,:,counter] = hist2d

        counter = counter+1

    hist_average_num = np.average(hist_array_num, axis = 2)
    hist_std_num = np.std(hist_array_num, axis = 2)
    hist_std_norm_num = hist_std_num / hist_average_num
    hist_std_norm_num = np.ma.masked_invalid(hist_std_norm_num)

    hist_average_mass = np.average(hist_array_mass, axis = 2)
    hist_std_mass = np.std(hist_array_mass, axis = 2)
    hist_std_norm_mass = hist_std_mass / hist_average_mass
    hist_std_norm_mass = np.ma.masked_invalid(hist_std_norm_mass)

    np.savetxt(f1, x_axis.centers())
    np.savetxt(f2, y_axis.centers())
    np.savetxt(f3, hist_average_num)
    np.savetxt(f4, hist_std_num)
    np.savetxt(f5, hist_average_mass)
    np.savetxt(f6, hist_std_norm_mass)

for case in ["10K_wei+1", "10K_wei-1", "10K_wei-4" ]:
    for hour in range(12, 13):
        files = []
        for i_loop in range (0, config.i_loop_max):
            filename_in = "urban_plume_wc_%s_0%03d_000000%02d.nc" % (case, (i_loop+1), hour)
            files.append(filename_in)
        f1 = "data/2d_scrit_%s_%02d_x_values.txt" % (case, hour)
        f2 = "data/2d_scrit_%s_%02d_y_values.txt" % (case, hour)
        f3 = "data/2d_scrit_%s_%02d_hist_average_num.txt" % (case, hour)
        f4 = "data/2d_scrit_%s_%02d_hist_std_norm_num.txt" % (case, hour)
        f5 = "data/2d_scrit_%s_%02d_hist_average_mass.txt" % (case, hour)
        f6 = "data/2d_scrit_%s_%02d_hist_std_norm_mass.txt" % (case, hour)
        print f1
        make_plot(files, f1, f2, f3, f4, f5, f6)




    



