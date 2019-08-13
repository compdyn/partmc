#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io, numpy

#axes2 = axes.twinx()

ncf = scipy.io.netcdf_file("out/out_mono_base/urban_plume_aq_chem_mono_process.nc")
time = ncf.variables["time"].data / 60
temp = ncf.variables["temp"].data
pressure = ncf.variables["pressure"].data

R_d = 287.058
conv_fac = pressure / R_d / temp
print(conv_fac[0])
print(ncf.variables["tot_num_conc"].data[0])

num_conc1 = ncf.variables["tot_num_conc"].data / conv_fac
dry_mass_conc1 = ncf.variables["tot_dry_mass_conc"].data / conv_fac
tot_mass_conc1 = ncf.variables["tot_mass_conc"].data / conv_fac
temp1 = temp

ncf = scipy.io.netcdf_file("out/out_mono_nh3/urban_plume_aq_chem_mono_process.nc")
time = ncf.variables["time"].data / 60
temp = ncf.variables["temp"].data
pressure = ncf.variables["pressure"].data

R_d = 287.058
conv_fac = pressure / R_d / temp
print(conv_fac[0])
print(ncf.variables["tot_num_conc"].data[0])

num_conc2 = ncf.variables["tot_num_conc"].data / conv_fac
dry_mass_conc2 = ncf.variables["tot_dry_mass_conc"].data / conv_fac
tot_mass_conc2 = ncf.variables["tot_mass_conc"].data / conv_fac
temp2 = temp

(figure, axes) = mpl_helper.make_fig(right_margin=2)
#axes.plot(time, num_conc, "b-")
#axes2.plot(time, dry_mass_conc, "r-")
axes.plot(time, dry_mass_conc1, "g-", label="base")
axes.plot(time, dry_mass_conc2, "r--", label="nh3")
axes.set_xlabel(r"time / min")
#axes.set_ylim([0,0.012])
#axes.set_ylabel(r"num. mix. ratio / $\rm kg^{-1}$")
axes.set_ylabel(r"dry mass mix. ratio / $\rm kg \, kg^{-1}$")
axes.grid(True)
axes.legend(loc=(1.05,0))
figure.savefig("figs/urban_plume_aq_chem_mono_total_dry.pdf")

(figure, axes) = mpl_helper.make_fig(right_margin=2)
axes.plot(time, tot_mass_conc1-dry_mass_conc1, "g-", label="base")
axes.plot(time, tot_mass_conc2-dry_mass_conc2, "r--", label="nh3")
axes.set_xlabel(r"time / min")
axes.set_ylabel(r"LWC / $\rm kg \, kg^{-1}$")
axes.grid(True)
axes.legend(loc=(1.05,0))
figure.savefig("figs/urban_plume_aq_chem_mono_total_h2o.pdf")

(figure, axes) = mpl_helper.make_fig(right_margin=2)
axes.plot(time, temp1, "g-", label="base")
axes.plot(time, temp2, "r--", label="nh3")
axes.set_xlabel(r"time / min")
axes.set_ylabel(r"temperature / K")
axes.grid(True)
axes.legend(loc=(1.05,0))
figure.savefig("figs/urban_plume_aq_chem_mono_temp.pdf")
