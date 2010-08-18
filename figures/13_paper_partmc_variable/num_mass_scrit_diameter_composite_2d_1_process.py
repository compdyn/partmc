#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
sys.path.append("../../tool")
import partmc
import config

i_loop_max = config.i_loop_max

def make_plot(hour, f1, f2, f3, f4, f5, f6):
    x_axis = partmc.log_grid(min=1e-9,max=1e-5,n_bin=70)
    y_axis = partmc.log_grid(min=1e-3,max=1e2,n_bin=50)
    x_centers = x_axis.centers()
    y_centers = y_axis.centers()

    hist_array_num = np.zeros([len(x_centers), len(y_centers), config.i_weighting_schemes, config.i_loop_max])
    hist_array_mass = np.zeros([len(x_centers), len(y_centers), config.i_weighting_schemes, config.i_loop_max])
    hist_average_num = np.zeros([len(x_centers), len(y_centers), config.i_weighting_schemes])
    hist_average_mass = np.zeros([len(x_centers), len(y_centers), config.i_weighting_schemes])
    hist_var_num = np.zeros([len(x_centers), len(y_centers), config.i_weighting_schemes])
    hist_var_mass = np.zeros([len(x_centers), len(y_centers), config.i_weighting_schemes])
    weighting_factor_num = np.zeros([len(x_centers), len(y_centers),config.i_weighting_schemes])
    weighting_factor_mass = np.zeros([len(x_centers), len(y_centers),config.i_weighting_schemes])

    for (counter_weighting, counter) in enumerate(["1K_wei+1", "1K_flat", "1K_wei-1", "1K_wei-2", "1K_wei-3", "1K_wei-4"]):
        print "I'm doing ", counter
        files = []
        for i_loop in range(0,config.i_loop_max):
            filename_in = config.netcdf_dir+"/urban_plume_wc_%s_0%03d_000000%02d.nc" % (counter,i_loop+1,hour)
            files.append(filename_in)

        for (counter_i_loop, file) in enumerate(files):
            print "file ", file
            ncf = Scientific.IO.NetCDF.NetCDFFile(dir_name+file)
            particles = partmc.aero_particle_array_t(ncf)
            env_state = partmc.env_state_t(ncf)
            ncf.close()

            dry_diameters = particles.dry_diameters()
            s_crit = (particles.critical_rel_humids(env_state) - 1)*100

            hist2d = partmc.histogram_2d(dry_diameters, s_crit, x_axis, y_axis, weights = 1/particles.comp_vols)
            hist_array_num[:,:,counter_weighting, counter_i_loop] = hist2d
            hist2d = partmc.histogram_2d(dry_diameters, s_crit, x_axis, y_axis, weights = particles.masses(include="[BC]")/particles.comp_vols)
            hist_array_mass[:,:,counter_weighting, counter_i_loop] = hist2d

    hist_array_num =  np.ma.masked_less_equal(hist_array_num,0)
    hist_array_mass =  np.ma.masked_less_equal(hist_array_mass,0)
    hist_average_num = np.average(hist_array_num, axis = 3)
    hist_average_mass = np.average(hist_array_mass, axis = 3)
    hist_var_num = np.var(hist_array_num, axis = 3)
    hist_var_mass = np.var(hist_array_mass, axis = 3)

    print "Calculated average and variance", counter

    weighting_factor_num = 1 / hist_var_num
    weighting_factor_mass = 1 / hist_var_mass
    weighting_factor_num_sum = np.sum(weighting_factor_num, axis = 2)
    weighting_factor_mass_sum = np.sum(weighting_factor_mass, axis = 2)

    hist_composite_num = np.zeros([len(x_centers), len(y_centers)])
    hist_composite_mass = np.zeros([len(x_centers), len(y_centers)])

    for i in range(0,config.i_weighting_schemes):
        increment = weighting_factor_num[:,:,i] / weighting_factor_num_sum * hist_average_num[:,:,i]
        increment = increment.filled(0)
        hist_composite_num += increment
    hist_composite_num = np.ma.masked_less_equal(hist_composite_num,0)

    for i in range(0,config.i_weighting_schemes):
        increment = weighting_factor_mass[:,:,i] / weighting_factor_mass_sum * hist_average_mass[:,:,i]
        increment = increment.filled(0)
        hist_composite_mass += increment
    hist_composite_mass = np.ma.masked_less_equal(hist_composite_mass,0)

np.savetxt(f1, x_axis.centers())
np.savetxt(f2, y_axis.centers())
np.savetxt(f3, hist_composite_num)
np.savetxt(f4, hist_composite_mass)

for hour in range(12, 13):
    f1 = "data/2d_scrit_compo_%02d_x_values.txt" % hour
    f2 = "data/2d_scrit_compo_%02d_y_values.txt" % hour
    f3 = "data/2d_scrit_compo_%02d_average_num.txt" % hour
    f4 = "data/2d_scrit_compo_%02d_average_mass.txt" % case, hour

    print f1
    make_plot(hour, f1, f2, f3, f4)




    



