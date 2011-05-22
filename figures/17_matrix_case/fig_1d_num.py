#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_dir, in_filename, out_filename, out_data_name):
    print in_filename
    ncf = scipy.io.netcdf.netcdf_file(in_dir+in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    x_axis = partmc.log_grid(min=1e-9,max=1e-5,n_bin=100)
    x_centers = x_axis.centers() 

    diameters = particles.diameters()

    pure_bc = ((particles.masses(include = ["BC"]) > 0) & (particles.masses(include = ["SO4"]) == 0))
    pure_so4 =  ((particles.masses(include = ["SO4"]) > 0) & (particles.masses(include = ["BC"]) == 0))
    with_bc = (particles.masses(include = ["BC"]) > 0)
    with_so4 = (particles.masses(include = ["SO4"]) > 0)
    mixed_bc_so4 = ((particles.masses(include = ["SO4"]) > 0) &  (particles.masses(include = ["BC"]) > 0))

    hist = partmc.histogram_1d(diameters, x_axis, weights = 1 / particles.comp_vols) / 1e6
    hist_bc = partmc.histogram_1d(diameters[pure_bc], x_axis, weights = 1 / particles.comp_vols[pure_bc]) /1e6
    hist_so4 = partmc.histogram_1d(diameters[pure_so4], x_axis, weights = 1 / particles.comp_vols[pure_so4]) /1e6
    hist_mixed = partmc.histogram_1d(diameters[mixed_bc_so4], x_axis, weights = 1 / particles.comp_vols[mixed_bc_so4]) / 1e6

    plt.clf()
    plt.loglog(x_centers*1e6, hist,  'r-', label = 'total')
    plt.loglog(x_centers*1e6, hist_bc,  'k-', label = 'pure bc')
    plt.loglog(x_centers*1e6, hist_so4,  'b-', label = 'pure so4')
    plt.loglog(x_centers*1e6, hist_mixed,  'g-', label = 'mixed so4 and bc')
    plt.axis([1e-3, 2e-0, 1e-1, 1e4])
    plt.xlabel("dry diameter / micrometer")
    plt.ylabel("number density / cm^{-3}")
    plt.legend(loc = "upper left")
    plt.grid(True)
    fig = plt.gcf()
    fig.savefig(out_filename)
    np.savetxt("diameter_values.txt", x_centers*1e6)
    np.savetxt(out_data_name+"_total_acc_bc1.txt", hist)
    np.savetxt(out_data_name+"_bc_acc_bc1.txt", hist_bc)
    np.savetxt(out_data_name+"_so4_acc_bc1.txt", hist_so4)
    np.savetxt(out_data_name+"_mixed_acc_bc1.txt", hist_mixed)

dir_name = "../../local_scenarios/matrix_case/out/"

filename_in = "brownian_part_0001_00000001.nc"
filename_out = "figs/1d_num_001_acc_bc1.pdf"
out_data_name = "data_001"
make_plot(dir_name, filename_in, filename_out, out_data_name)

filename_in = "brownian_part_0001_00000073.nc"
filename_out = "figs/1d_num_073_acc_bc1.pdf"
out_data_name = "data_073"
make_plot(dir_name, filename_in, filename_out, out_data_name)

filename_in = "brownian_part_0001_00000145.nc"
filename_out = "figs/1d_num_145_acc_bc1.pdf"
out_data_name = "data_145"
make_plot(dir_name, filename_in, filename_out, out_data_name)
