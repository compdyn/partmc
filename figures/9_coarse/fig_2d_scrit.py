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

def make_plot(hour_counter, case_counter, dir_name,in_files,out_filename1, out_filename2):
    x_axis = partmc.log_grid(min=1e-9,max=1e-5,n_bin=70)
    y_axis = partmc.log_grid(min=1e-3,max=1e2,n_bin=50)
    x_centers = x_axis.centers()
    y_centers = y_axis.centers()
    i_counter = 0
    hist_array = np.zeros([len(x_centers), len(y_centers), config.i_loop_max])
    hist_average = np.zeros([len(x_centers), len(y_centers)])
    hist_std = np.zeros([len(x_centers), len(y_centers)])
    hist_std_norm = np.zeros([len(x_centers), len(y_centers)])
    ccn_array = np.zeros([3,config.i_loop_max])

    for file in in_files:
        ncf = Scientific.IO.NetCDF.NetCDFFile(dir_name+file)
        particles = partmc.aero_particle_array_t(ncf)
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        dry_diameters = particles.dry_diameters()
        total_number = len(particles.ids)
        s_crit = (particles.critical_rel_humids(env_state) - 1)*100
        hist2d = partmc.histogram_2d(dry_diameters, s_crit, x_axis, y_axis, weights = 1/particles.comp_vols)
        hist_array[:,:,i_counter] = hist2d

        activated_1 = (s_crit < config.s_crit_1)
        number_act_1 = sum(1/particles.comp_vols[activated_1])

        activated_2 = (s_crit < config.s_crit_2)
        number_act_2 = sum(1/particles.comp_vols[activated_2])

        activated_3 = (s_crit < config.s_crit_3)
        number_act_3 = sum(1/particles.comp_vols[activated_3])

        ccn_array[0,i_counter]= number_act_1/total_number
        ccn_array[1,i_counter]= number_act_2/total_number
        ccn_array[2,i_counter]= number_act_3/total_number

        i_counter += 1

    hist_average = np.average(hist_array, axis = 2)
    hist_std = np.std(hist_array, axis = 2)
    hist_std_norm = hist_std/hist_average
    hist_std_norm = np.nan_to_num(hist_std_norm)
    
    print 'counters 3 ', hour_counter, case_counter

    ccn = np.average(ccn_array, axis = 1)
    ccn_average[hour_counter,case_counter,:] = ccn

    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist_average.transpose(), norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("log")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("critical supersaturation (%)")
    plt.clim(1e8, 1e12)
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    plt.grid()
    fig = plt.gcf()
    fig.savefig(out_filename1)

    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist_std_norm.transpose(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("log")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("critical supersaturation (%)")
    cbar = plt.colorbar()
    plt.clim(0, 3)
    cbar.set_label("std/avg")
    plt.grid()
    fig = plt.gcf()
    fig.savefig(out_filename2)

dir_name = "../../scenarios/5_weighted/out_10loop/"
hour_counter = 0
ccn_average = np.zeros([25,10,3])

for hour in range(1,26):
    print "hour = ", hour
    case_counter = 0
    print 'counters 1 ', hour_counter, case_counter
    for counter in ["10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4", "100K_flat", "100K_wei-1", "100K_wei-2", "100K_wei-3", "100K_wei-4"]:
        print 'counter ', counter
        files = []
        for i_loop in range(0,config.i_loop_max):
            filename_in = "urban_plume_wc_%s_00%02d_000000%02d.nc" % (counter,i_loop+1,hour)
            files.append(filename_in)
        filename_out1 = "figs/2d_scrit_%s_%02d.pdf" % (counter, hour)
        filename_out2 = "figs/2d_scrit_std_%s_%02d.pdf" % (counter, hour)
        make_plot(hour_counter, case_counter, dir_name, files, filename_out1, filename_out2)
        print 'counters 2', hour_counter, case_counter
        case_counter = case_counter + 1
    hour_counter = hour_counter + 1
np.savetxt("data/ccn_average_ss1.txt", ccn_average[:,:,0])
np.savetxt("data/ccn_average_ss2.txt", ccn_average[:,:,1])
np.savetxt("data/ccn_average_ss3.txt", ccn_average[:,:,2])

#    for counter in ["10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4", "100K_flat", "100K_wei-1", "100K_wei-2", "100K_wei-3", "100K_wei-4"]:
#    for counter in ["10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4"]:


