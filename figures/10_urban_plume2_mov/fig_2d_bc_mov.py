#!/usr/bin/env python

import config_matplotlib
import scipy.io
import sys
import numpy as np
import matplotlib
#matplotlib.use("PDF")
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

#matplotlib.rc('text', usetex = True)
#matplotlib.rc('font', size = 20, family = "serif",
#              serif = ["Computer Modern Roman"])
#matplotlib.rc('xtick.major', pad = 5)
#matplotlib.rc('ytick.major', pad = 5)
#matplotlib.rc('xtick', labelsize = 20)
#matplotlib.rc('ytick', labelsize = 20)
#

def make_plot(in_filename,out_filename):
    ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    bc = particles.masses(include = ["BC"])
    dry_mass = particles.masses(exclude = ["H2O"])
    bc_frac = bc / dry_mass

    wet_diameters = particles.diameters()
    dry_diameters = particles.dry_diameters() * 1e6

    x_axis = partmc.log_grid(min=1e-2,max=1e0,n_bin=90)
    y_axis = partmc.linear_grid(min=0,max=0.8,n_bin=40)

    (figure, axes, cbar_axes) = config_matplotlib.make_fig(colorbar=True,
                                                           right_margin=1,
                                                           top_margin=0.3)

    axes.grid(True)
    axes.grid(True, which = 'minor')
    axes.minorticks_on()
    axes.set_xscale('log')

    axes.set_xbound(x_axis.min, x_axis.max)
    axes.set_ybound(y_axis.min, y_axis.max)

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8
    #xaxis.set_major_formatter(matplotlib.ticker.LogFormatter())
    #yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    #yaxis.set_minor_locator(matplotlib.ticker.MaxNLocator(8))

    axes.set_xlabel(r"dry diameter $D\ (\rm\mu m)$")
    axes.set_ylabel(r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}$")

    hist2d = partmc.histogram_2d(dry_diameters, bc_frac, x_axis, y_axis, weights = 1/particles.comp_vols)
#    plt.clf()

    axes.set_xbound(x_axis.min, x_axis.max)
    axes.set_ybound(y_axis.min, y_axis.max)

    p=axes.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(vmin=1e8, vmax=1e11), linewidths = 0.1)
    title = partmc.time_of_day_string(env_state)
    axes.set_xbound(x_axis.min, x_axis.max)
    axes.set_ybound(y_axis.min, y_axis.max)
    axes.set_title(title)
    figure.colorbar(p, cax = cbar_axes, format = matplotlib.ticker.LogFormatterMathtext())
    cbar_axes.set_ylabel(r"number conc. $(\rm m^{-3})$")
    #cbar_axes.set_ylim([1e8, 1e11])
    #plt.title(title)
    axes.set_xbound(x_axis.min, x_axis.max)
    axes.set_ybound(y_axis.min, y_axis.max)
    fig = plt.gcf()
    fig.savefig(out_filename)

for counter in range(1, 1442):
    print "counter = ",  counter
    
    filename_in1 = "../../scenarios/2_urban_plume2/out/urban_plume_wc_0001_0000%04d.nc" % counter
    filename_out1 = "figs/2d_bc_%04d.pdf" % (counter-1)
    print filename_in1
    print filename_out1

    make_plot(filename_in1, filename_out1)
#    make_plot(filename_in2, filename_out2, titel)
#    make_plot(filename_in3, filename_out3, titel)
#    make_plot(filename_in4, filename_out4, titel)

