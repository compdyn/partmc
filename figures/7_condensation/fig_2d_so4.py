#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

def make_plot(in_filename,out_filename,title):
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

    so4 = particles.mass(include = ["SO4"])/particles.aero_data.molec_weight[0]
    nh4 =  particles.mass(include = ["NH4"])/particles.aero_data.molec_weight[3]
    no3 =  particles.mass(include = ["NO3"])/particles.aero_data.molec_weight[1]

    dry_mass = particles.mass(exclude = ["H2O"])
    so4_frac = so4 / dry_mass
    ion_ratio = (2*so4 + no3) / nh4

    is_neutral = (ion_ratio < 2)
    print 'neutral ', sum(is_neutral), ion_ratio[is_neutral]

    dry_diameter = particles.dry_diameter()

    x_axis = pmc_data_nc.pmc_log_axis(min=1e-8,max=1e-6,n_bin=70)
    y_axis = pmc_data_nc.pmc_linear_axis(min=0,max=1.0,n_bin=50)

    hist2d = pmc_data_nc.histogram_2d(dry_diameter, so4_frac, x_axis, y_axis, weights = 1/particles.comp_vol)

    plt.clf()
    plt.semilogx(dry_diameter, ion_ratio, 'rx')
    fig = plt.gcf()
    fig.savefig('figs/t.pdf')

    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("SO4 mass fraction")
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for hour in range(2, 3):
    print "hour = ", hour
    
    filename_in1 = "../../urban_plume2/out_no_nh3/urban_plume_nc_0001_000000%02d.nc" % hour
    filename_out1 = "figs/2d_so4_no_nh3_nc_%02d.pdf" % (hour-1)
    titel = "%02d hours" % (hour-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, titel)


