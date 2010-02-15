#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
const = partmc.constants_t("../../src/constants.f90")

def make_plot(in_filename,out_filename,title):
    print in_filename
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename)
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    dry_diameters = particles.dry_diameters()
    s_crit = (particles.critical_rel_humids(env_state, const) - 1)*100
    x_axis = partmc.log_grid(min=1e-8,max=1e-6,n_bin=70)
    y_axis = partmc.log_grid(min=1e-3,max=1e2,n_bin=50)

    hist2d = partmc.histogram_2d(dry_diameters, s_crit, x_axis, y_axis, weights = 1/particles.comp_vols)
    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("log")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("critical supersaturation (%)")
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for hour in range(1, 49):
    print "hour = ", hour

    filename_in1 = "../../scenarios/3_condense/start/urban_plume_wc_0001_000000%02d.nc" % hour
    filename_in2 = "../../scenarios/3_condense/start/urban_plume_comp_wc_0001_000000%02d.nc" % hour
    filename_in3 = "../../scenarios/3_condense/start/urban_plume_size_wc_0001_000000%02d.nc" % hour
    filename_in4 = "../../scenarios/3_condense/start/urban_plume_both_wc_0001_000000%02d.nc" % hour
    filename_out1 = "figs/2d_scrit_ref_%02d.pdf" % (hour-1)
    filename_out2 = "figs/2d_scrit_comp_%02d.pdf" % (hour-1)
    filename_out3 = "figs/2d_scrit_size_%02d.pdf" % (hour-1)
    filename_out4 = "figs/2d_scrit_both_%02d.pdf" % (hour-1)
    titel = "%02d hours" % (hour-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, titel)
    make_plot(filename_in2, filename_out2, titel)
    make_plot(filename_in3, filename_out3, titel)
    make_plot(filename_in4, filename_out4, titel)
 
#make_plot("../../scenarios/3_condense/start_up2/urban_plume_wc_0001_00000031.nc","figs/2d_scrit_30.pdf","30 hours")
#make_plot("../../scenarios/3_condense/start_up2/urban_plume_wc_0001_00000037.nc","figs/2d_scrit_36.pdf","36 hours")
#make_plot("../../scenarios/3_condense/start_up2/urban_plume_wc_0001_00000043.nc","figs/2d_scrit_42.pdf","42 hours")
#make_plot("../../scenarios/3_condense/start_up2/urban_plume_wc_0001_00000049.nc","figs/2d_scrit_48.pdf","48 hours")

#make_plot("../../scenarios/3_condense/start_up2_comp/urban_plume_wc_0001_00000031.nc","figs/2d_scrit_comp_30.pdf","30 hours")
#make_plot("../../scenarios/3_condense/start_up2_comp/urban_plume_wc_0001_00000037.nc","figs/2d_scrit_comp_36.pdf","36 hours")
#make_plot("../../scenarios/3_condense/start_up2_comp/urban_plume_wc_0001_00000043.nc","figs/2d_scrit_comp_42.pdf","42 hours")
#make_plot("../../scenarios/3_condense/start_up2_comp/urban_plume_wc_0001_00000049.nc","figs/2d_scrit_comp_48.pdf","48 hours")

#make_plot("../../scenarios/3_condense/start_up2_size/urban_plume_wc_0001_00000031.nc","figs/2d_scrit_size_30.pdf","30 hours")
#make_plot("../../scenarios/3_condense/start_up2_size/urban_plume_wc_0001_00000037.nc","figs/2d_scrit_size_36.pdf","36 hours")
#make_plot("../../scenarios/3_condense/start_up2_size/urban_plume_wc_0001_00000043.nc","figs/2d_scrit_size_42.pdf","42 hours")
#make_plot("../../scenarios/3_condense/start_up2_size/urban_plume_wc_0001_00000049.nc","figs/2d_scrit_size_48.pdf","48 hours")

#make_plot("../../scenarios/3_condense/start_up2_both/urban_plume_wc_0001_00000031.nc","figs/2d_scrit_both_30.pdf","30 hours")
#make_plot("../../scenarios/3_condense/start_up2_both/urban_plume_wc_0001_00000037.nc","figs/2d_scrit_both_36.pdf","36 hours")
#make_plot("../../scenarios/3_condense/start_up2_both/urban_plume_wc_0001_00000043.nc","figs/2d_scrit_both_42.pdf","42 hours")
#make_plot("../../scenarios/3_condense/start_up2_both/urban_plume_wc_0001_00000049.nc","figs/2d_scrit_both_48.pdf","48 hours")
