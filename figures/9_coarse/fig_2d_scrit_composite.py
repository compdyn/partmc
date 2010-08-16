#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import config

def make_plot(dir_name, hour, out_filename1):
    x_axis = partmc.log_grid(min=1e-9,max=1e-5,n_bin=70)
    y_axis = partmc.log_grid(min=1e-3,max=1e2,n_bin=50)
    x_centers = x_axis.centers()
    y_centers = y_axis.centers()
   
    hist_array = np.zeros([len(x_centers), len(y_centers), config.i_weighting_schemes, config.i_loop_max])
    hist_average = np.zeros([len(x_centers), len(y_centers), config.i_weighting_schemes])
    hist_var = np.zeros([len(x_centers), len(y_centers), config.i_weighting_schemes])
    weighting_factor = np.zeros([len(x_centers), len(y_centers),config.i_weighting_schemes])
    
    for (counter_weighting, counter) in enumerate(["10K_wei+1", "10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4"]):
        print "I'm doing ", counter
        files = []
        for i_loop in range(0,config.i_loop_max):
            filename_in = "urban_plume_wc_%s_0%03d_000000%02d.nc" % (counter,i_loop+1,hour)
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
            hist_array[:,:,counter_weighting, counter_i_loop] = hist2d
           
    hist_array =  np.ma.masked_less_equal(hist_array,0)
    hist_average = np.average(hist_array, axis = 3)
    hist_var = np.var(hist_array, axis = 3)
    print "Calculated average and variance", counter

    weighting_factor = 1 / hist_var
    weighting_factor_sum = np.sum(weighting_factor, axis = 2)

    hist_composite = np.zeros([len(x_centers), len(y_centers)])
    for i in range(0,config.i_weighting_schemes):
        increment = weighting_factor[:,:,i] / weighting_factor_sum * hist_average[:,:,i]
        increment = increment.filled(0)
        hist_composite += increment
    hist_composite = np.ma.masked_less_equal(hist_composite,0)        

    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist_composite.transpose(), norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("log")
    plt.grid()
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("S_crit in %")
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    fig = plt.gcf()
    fig.savefig(out_filename1)


dir_name = "../../scenarios/5_weighted/out/"
for hour in range(12,13):
    print "hour = ", hour
    filename_out1 = "figs/2d_scrit_composite_10K_%02d.pdf" % hour
    make_plot(dir_name, hour, filename_out1)


