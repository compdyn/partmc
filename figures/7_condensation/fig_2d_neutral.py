#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
matplotlib.use('Agg')
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
    bc =  particles.mass(include = ["BC"])/particles.aero_data.molec_weight[18]
    oc =  particles.mass(include = ["OC"])/particles.aero_data.molec_weight[17]

    print 'min nh4 ', min(particles.mass(include = ["NH4"])), max(nh4), min(no3), max(no3)

    ion_ratio = (2*so4 + no3) / nh4

    is_neutral = (ion_ratio < 2)
    dry_diameter = particles.dry_diameter()

    x_axis = pmc_data_nc.pmc_log_axis(min=1e-8,max=1e-6,n_bin=70)
    y_axis = pmc_data_nc.pmc_linear_axis(min=0,max=30.0,n_bin=100)
    x_centers = x_axis.centers()

    bin_so4 = pmc_data_nc.histogram_1d(dry_diameter, x_axis, weights = so4)
    bin_nh4 = pmc_data_nc.histogram_1d(dry_diameter, x_axis, weights = nh4)
    bin_no3 = pmc_data_nc.histogram_1d(dry_diameter, x_axis, weights = no3)
    
    print 'bin_so4 ', bin_so4[40]
    print 'bin_nh4 ', bin_nh4[40]
    print 'bin_no3 ', bin_no3[40]

    bin_ratio = (2*bin_so4 + bin_no3)/ bin_nh4
    np.isnan(bin_ratio) # checks which elements in c are NaN (produces array with True and False)
    bin_ratio[np.isnan(bin_ratio)] = 0 # replaces NaN with 0. useful for plotting
    print 'bin_ratio ', bin_ratio[40]

    diameter_bins = x_axis.find(dry_diameter)
    print 'diameter_bins ', diameter_bins
    is_40 = (diameter_bins == 40)
#    for i in range(len(dry_diameter)):
#        if diameter_bins[i] == 40:
#            print 'particle info', so4[i], nh4[i], no3[i], ion_ratio[i]
    so4_40 = so4[is_40]
    nh4_40 = nh4[is_40]
    no3_40 = no3[is_40]
    bc_40 = bc[is_40]
    oc_40 = oc[is_40]

    ion_ratio_40 = ion_ratio[is_40]
#    data = [(so4_40[i],nh4_40[i], no3_40[i], ion_ratio_40[i]) for i in range(len(so4_40)
    data = zip(so4_40, nh4_40, no3_40, bc_40, oc_40, ion_ratio_40)
    data.sort(key = lambda x: x[5])
    for (so,nh,no,bc,oc,ir) in data:
        print so,nh,no,bc,oc,ir

    print 'sums ', sum(so4[is_40]), sum(nh4[is_40]), sum(no3[is_40]), (2*sum(so4[is_40])+ sum(no3[is_40])) / sum(nh4[is_40])
    print 'sums/number ',  sum(so4[is_40])/len(so4_40), sum(nh4[is_40])/len(nh4_40), sum(no3[is_40])/len(no3_40)
    
    
    hist2d = pmc_data_nc.histogram_2d(dry_diameter, ion_ratio, x_axis, y_axis, weights = 1/particles.comp_vol)

    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    plt.semilogx(x_centers, bin_ratio, 'w-', linewidth = 3)
    plt.semilogx(x_centers, bin_ratio, 'k-', linewidth = 1)
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("ion ratio")
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for hour in range(13, 14):
    print "hour = ", hour
    
    filename_in1 = "../../urban_plume2/out_no_nh3/urban_plume_wc_0001_000000%02d.nc" % hour
    filename_out1 = "figs/2d_neutral_no_nh3_wc_%02d.png" % (hour-1)
    titel = "%02d hours" % (hour-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, titel)


