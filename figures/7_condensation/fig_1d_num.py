#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

ccn_cn_ratio = np.zeros([4,49])
def make_plot(in_dir, in_filename1, in_filename2, in_filename3, out_filename, title, ccn_cn_i, ccn_cn_j):
    print in_filename1, in_filename2, in_filename3
    ncf = scipy.io.netcdf.netcdf_file(in_dir+in_filename1, 'r')
    particles1 = partmc.aero_particle_array_t(ncf)
    ncf.close()
    ncf = scipy.io.netcdf.netcdf_file(in_dir+in_filename2, 'r')
    particles2 = partmc.aero_particle_array_t(ncf)
    ncf.close()
    ncf = scipy.io.netcdf.netcdf_file(in_dir+in_filename3, 'r')
    particles3 = partmc.aero_particle_array_t(ncf)
    ncf.close()

    x_axis = partmc.log_grid(min=1e-10,max=1e-4,n_bin=30)
    x_centers = x_axis.centers() 

    wet_diameters1 = particles1.diameters()
    wet_diameters2 = particles2.diameters()
    wet_diameters3 = particles3.diameters()

    hist1 = partmc.histogram_1d(wet_diameters1, x_axis, weights = 1 / particles1.comp_vols)
    hist2 = partmc.histogram_1d(wet_diameters2, x_axis, weights = 1 / particles2.comp_vols)
    hist3 = partmc.histogram_1d(wet_diameters3, x_axis, weights = 1 / particles3.comp_vols)
 
    is_activated = (wet_diameters3 > 2e-6)
    sum_tot = sum(1/particles3.comp_vols) * 1e-6
    num_act = sum(1/particles3.comp_vols[is_activated]) * 1e-6
    print title, num_act, sum_tot, num_act/sum_tot * 100

    ccn_cn_ratio[ccn_cn_i, ccn_cn_j] =  num_act/sum_tot

    plt.clf()
    plt.semilogx(x_axis.centers(), hist1, label = '0 min')
    plt.semilogx(x_axis.centers(), hist2, label = '2 mins')
    plt.semilogx(x_axis.centers(), hist3, label = '10 mins') 
    plt.legend(loc = 'upper left')
    plt.xlabel("wet diameter (m)")
    plt.ylabel("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for counter in range(1,41):
    print "counter = ", counter

    dir_name = "../../scenarios/3_condense/out/"

    filename_in1 = "cond_%02d_ref_0001_00000001.nc" % counter
    filename_in2 = "cond_%02d_ref_0001_00000121.nc" % counter
    filename_in3 = "cond_%02d_ref_0001_00000601.nc" % counter

    filename_out = "figs/1d_wc_num_%02d_ref.pdf" % (counter - 1)
    title = "%02d hours" % (counter - 1)

    print dir_name, title
    print filename_in1, filename_in2, filename_in3
    print filename_out

    make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out, title, 0, counter-1)
    filename_in1 = "cond_%02d_comp_0001_00000001.nc" % counter
    filename_in2 = "cond_%02d_comp_0001_00000121.nc" % counter
    filename_in3 = "cond_%02d_comp_0001_00000601.nc" % counter

    filename_out = "figs/1d_wc_num_%02d_comp.pdf" % (counter - 1)
    title = "%02d hours" % (counter - 1)

    print dir_name, title
    print filename_in1, filename_in2, filename_in3
    print filename_out

    make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out, title, 1, counter-1)

    filename_in1 = "cond_%02d_size_0001_00000001.nc" % counter
    filename_in2 = "cond_%02d_size_0001_00000121.nc" % counter
    filename_in3 = "cond_%02d_size_0001_00000601.nc" % counter

    filename_out = "figs/1d_wc_num_%02d_size.pdf" % (counter - 1)
    title = "%02d hours" % (counter - 1)

    print dir_name, title
    print filename_in1, filename_in2, filename_in3
    print filename_out

    make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out, title, 2, counter-1)
    
    filename_in1 = "cond_%02d_both_0001_00000001.nc" % counter
    filename_in2 = "cond_%02d_both_0001_00000121.nc" % counter
    filename_in3 = "cond_%02d_both_0001_00000601.nc" % counter

    filename_out = "figs/1d_wc_num_%02d_both.pdf" % (counter - 1)
    title = "%02d hours" % (counter - 1)

    print dir_name, title
    print filename_in1, filename_in2, filename_in3
    print filename_out

    make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out, title, 3, counter-1)

np.savetxt("data/ccn_cn_ratio_wc.txt", ccn_cn_ratio)
