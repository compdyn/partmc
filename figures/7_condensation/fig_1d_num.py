#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

ccn_cn_ratio = np.zeros([4,49])
def make_plot(in_dir, in_filename1, in_filename2, in_filename3, out_filename, title, ccn_cn_i, ccn_cn_j):
    print in_filename1, in_filename2, in_filename3
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename1)
    particles1 = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename2)
    particles2 = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename3)
    particles3 = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

    x_axis = pmc_data_nc.pmc_log_axis(min=1e-10,max=1e-4,n_bin=30)
    x_centers = x_axis.centers() 

    wet_diameter1 = particles1.diameter()
    wet_diameter2 = particles2.diameter()
    wet_diameter3 = particles3.diameter()

    hist1 = pmc_data_nc.histogram_1d(wet_diameter1, x_axis, weights = 1 / particles1.comp_vol)
    hist2 = pmc_data_nc.histogram_1d(wet_diameter2, x_axis, weights = 1 / particles2.comp_vol)
    hist3 = pmc_data_nc.histogram_1d(wet_diameter3, x_axis, weights = 1 / particles3.comp_vol)
 
    is_activated = (wet_diameter3 > 2e-6)
    sum_tot = sum(1/particles3.comp_vol) * 1e-6
    num_act = sum(1/particles3.comp_vol[is_activated]) * 1e-6
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

    dir_name = "../../new_cond/out/"

    filename_in1 = "cond_wc_%02d_ref_0001_00000001.nc" % counter
    filename_in2 = "cond_wc_%02d_ref_0001_00000121.nc" % counter
    filename_in3 = "cond_wc_%02d_ref_0001_00000601.nc" % counter

    filename_out = "figs/1d_wc_num_%02d_ref.pdf" % (counter - 1)
    title = "%02d hours" % (counter - 1)

    print dir_name, title
    print filename_in1, filename_in2, filename_in3
    print filename_out

    make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out, title, 0, counter-1)
    filename_in1 = "cond_wc_%02d_comp_0001_00000001.nc" % counter
    filename_in2 = "cond_wc_%02d_comp_0001_00000121.nc" % counter
    filename_in3 = "cond_wc_%02d_comp_0001_00000601.nc" % counter

    filename_out = "figs/1d_wc_num_%02d_comp.pdf" % (counter - 1)
    title = "%02d hours" % (counter - 1)

    print dir_name, title
    print filename_in1, filename_in2, filename_in3
    print filename_out

    make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out, title, 1, counter-1)

    filename_in1 = "cond_wc_%02d_size_0001_00000001.nc" % counter
    filename_in2 = "cond_wc_%02d_size_0001_00000121.nc" % counter
    filename_in3 = "cond_wc_%02d_size_0001_00000601.nc" % counter

    filename_out = "figs/1d_wc_num_%02d_size.pdf" % (counter - 1)
    title = "%02d hours" % (counter - 1)

    print dir_name, title
    print filename_in1, filename_in2, filename_in3
    print filename_out

    make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out, title, 2, counter-1)
    
    filename_in1 = "cond_wc_%02d_both_0001_00000001.nc" % counter
    filename_in2 = "cond_wc_%02d_both_0001_00000121.nc" % counter
    filename_in3 = "cond_wc_%02d_both_0001_00000601.nc" % counter

    filename_out = "figs/1d_wc_num_%02d_both.pdf" % (counter - 1)
    title = "%02d hours" % (counter - 1)

    print dir_name, title
    print filename_in1, filename_in2, filename_in3
    print filename_out

    make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out, title, 3, counter-1)

np.savetxt("data/ccn_cn_ratio_wc.txt", ccn_cn_ratio)
