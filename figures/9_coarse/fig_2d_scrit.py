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

def make_plot(dir_name,in_files,out_filename1, out_filename2):
    x_axis = partmc.log_grid(min=1e-9,max=1e-5,n_bin=70)
    y_axis = partmc.log_grid(min=1e-3,max=1e2,n_bin=50)
    x_centers = x_axis.centers()
    y_centers = y_axis.centers()
    counter = 0
    hist_array = np.zeros([len(x_centers), len(y_centers), config.i_loop_max])
    hist_average = np.zeros([len(x_centers), len(y_centers)])
    hist_std = np.zeros([len(x_centers), len(y_centers)])
    hist_std_norm = np.zeros([len(x_centers), len(y_centers)])

    for file in in_files:
        ncf = Scientific.IO.NetCDF.NetCDFFile(dir_name+file)
        particles = partmc.aero_particle_array_t(ncf)
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        dry_diameters = particles.dry_diameters()
        s_crit = (particles.critical_rel_humids(env_state) - 1)*100
        hist2d = partmc.histogram_2d(dry_diameters, s_crit, x_axis, y_axis, weights = 1/particles.comp_vols)
        hist_array[:,:,counter] = hist2d
        counter = counter + 1

    hist_average = np.average(hist_array, axis = 2)
    hist_std = np.std(hist_array, axis = 2)
    hist_std_norm = hist_std/hist_average
    hist_std_norm = np.ma.masked_invalid(hist_std_norm)

    print "min, max", hist_average.min(), hist_average.max()
  
#    dry_diameters_line = np.array([1e-9, 1e-5])
#    kappa_line1 = np.array([0.01, 0.01])
#    crit_ss_line1 = (partmc.critical_rel_humids(env_state,kappa_line1, dry_diameters_line)-1)*100.
#    kappa_line2 = np.array([0.1, 0.1])
#    crit_ss_line2 = (partmc.critical_rel_humids(env_state,kappa_line2, dry_diameters_line)-1)*100.
#    kappa_line3 = np.array([2,2])
#    crit_ss_line3 = (partmc.critical_rel_humids(env_state,kappa_line3, dry_diameters_line)-1)*100.
#    kappa_line4 = np.array([0.001,0.001])
#    crit_ss_line4 = (partmc.critical_rel_humids(env_state,kappa_line4, dry_diameters_line)-1)*100.
 
#    print 'line ', kappa_line1, dry_diameters_line, crit_ss_line1

    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist_average.transpose(), norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e11), linewidths = 0.1)
#    plt.plot(dry_diameters_line, crit_ss_line1, 'k-')
#    plt.plot(dry_diameters_line, crit_ss_line2, 'k-')
#    plt.plot(dry_diameters_line, crit_ss_line3, 'k-')
#    plt.plot(dry_diameters_line, crit_ss_line4, 'k-')
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("log")
    plt.grid()
    plt.axis([5e-9, 5e-6, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("S_crit in %")
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    fig = plt.gcf()
    fig.savefig(out_filename1)

    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist_std_norm.transpose(), norm = matplotlib.colors.LogNorm(vmin=1e-2, vmax = 10), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("log")
    plt.grid()
    plt.axis([5e-9, 5e-6, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("S_crit in %")
    cbar = plt.colorbar()
    cbar.set_label("CV")
    fig = plt.gcf()
    fig.savefig(out_filename2)

dir_name = "../../scenarios/5_weighted/out/"
for hour in range(12,13):
    print "hour = ", hour
#    for counter in ["1K_wei+1", "1K_flat", "1K_wei-1", "1K_wei-2", "1K_wei-3", "1K_wei-4", "10K_wei+1", "10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4"]:
    for counter in ["10K_wei+1", "10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4"]:
#    for counter in ["1K_flat"]:
        print 'counter ', counter
        files = []
        for i_loop in range(0,config.i_loop_max):
            filename_in = "urban_plume_wc_%s_0%03d_000000%02d.nc" % (counter,i_loop+1,hour)
            files.append(filename_in)
        filename_out1 = "figs/2d_scrit_%s_%02d.pdf" % (counter, hour)
        filename_out2 = "figs/2d_scrit_std_%s_%02d.pdf" % (counter, hour)
        make_plot(dir_name, files, filename_out1, filename_out2)


