#!/usr/bin/env python

import scipy.io
import sys
import math
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_filename,out_filename,title):
    ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    bc_volume = particles.volumes(include = ["BC"])
    bc  = particles.masses(include = ["BC"])
    dry_mass = particles.masses(exclude = ["H2O"])
    bc_frac = bc / dry_mass
    coat_frac = 1 - bc / dry_mass

    is_bc = (bc_frac > 0.05)
    dry_diameters = particles.dry_diameters()

    core_diameters = (6 / math.pi * bc_volume)**(1./3.)
    coating_thickness = (dry_diameters - core_diameters) / 2.
    ratio = coating_thickness/dry_diameters

    print ratio.max()
    x_axis = partmc.linear_grid(min=0,max=ratio.max(),n_bin=50)

    hist1d = partmc.histogram_1d(coating_thickness[is_bc]/dry_diameters[is_bc], x_axis, weights = 1/particles.comp_vols[is_bc])
    print hist1d
    plt.clf()
    a = plt.gca()
    a.set_xscale("linear")
    a.set_yscale("log")
    plt.plot(x_axis.centers(), hist1d)
    plt.axis([x_axis.min, x_axis.max, 1e8, 1e12])
    plt.grid(True)
    plt.xlabel("coating thickness / dry diameter")
    plt.ylabel("number concentration (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for hour in range(2, 26):
    minutes = hour 
    print "hour = ", hour
    filename_in1 = "../../local_scenarios/aging_comp/run_100K_60min/out/urban_plume2_wc_0001_0000%04d.nc" % minutes
    filename_out1 = "figs/coating_%02d.pdf" % (hour-1)
    titel = "%02d hours" % (hour-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, titel)
#    make_plot(filename_in2, filename_out2, titel)
#    make_plot(filename_in3, filename_out3, titel)
#    make_plot(filename_in4, filename_out4, titel)


