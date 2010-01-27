#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
const = pmc_data_nc.load_constants("../../src/constants.f90")

def make_plot(in_filename,out_filename,time,title):
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    env_state = pmc_data_nc.env_state_t(ncf)
    ncf.close()

    age = abs(particles.least_create_time / 3600. - time)
    dry_diameter = particles.dry_diameter()
    s_crit = (particles.critical_rh(env_state, const) - 1)*100

    x_axis = pmc_data_nc.pmc_log_axis(min=1e-8,max=1e-6,n_bin=140)
    y_axis = pmc_data_nc.pmc_linear_axis(min=0, max = 48, n_bin=96)

    vals2d = pmc_data_nc.multival_2d(dry_diameter, age, s_crit, x_axis, y_axis)

    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), vals2d.transpose(), norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("age (h)")
    cbar = plt.colorbar()
    cbar.set_label("S_crit (%)")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for hour in range(49, 50):
    print "hour = ", hour
    time = hour - 1
    filename_in1 = "../../urban_plume2/out_no_nh3/urban_plume_nc_0001_000000%02d.nc" % hour
    filename_out1 = "figs/2d_age_scrit_no_nh3_nc_%02d.pdf" % (hour-1)
    titel = "%02d hours" % (hour-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, time, titel)


