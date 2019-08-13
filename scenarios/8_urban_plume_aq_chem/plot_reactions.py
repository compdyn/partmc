#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io
import numpy as np

ncf = scipy.io.netcdf_file("out/out_mono_base/urban_plume_aq_chem_mono_process.nc")
time = ncf.variables["time"].data / 60

forw = []
back = []
conv = []
for i in range(1,602):
    ncf = scipy.io.netcdf_file(f"out/out_mono_base/urban_plume_aq_chem_mono_0001_{i:08d}.nc")
    forw.append(ncf.variables["aq_chem_rates_forward"].data.copy())
    back.append(ncf.variables["aq_chem_rates_backward"].data.copy())
    conv.append(ncf.variables["aq_chem_rates_conv_factors"].data.copy())
    ncf.close()


forw1 = np.array(forw) # time x reaction x particle
back1 = np.array(back)
conv1 = np.array(conv)

print(conv1[20,0:50,1])

ncf = scipy.io.netcdf_file("out/out_mono_nh3/urban_plume_aq_chem_mono_process.nc")
time = ncf.variables["time"].data / 60

forw = []
back = []
conv = []
for i in range(1,602):
    ncf = scipy.io.netcdf_file(f"out/out_mono_nh3/urban_plume_aq_chem_mono_0001_{i:08d}.nc")
    forw.append(ncf.variables["aq_chem_rates_forward"].data.copy())
    back.append(ncf.variables["aq_chem_rates_backward"].data.copy())
    conv.append(ncf.variables["aq_chem_rates_conv_factors"].data.copy())
    ncf.close()


forw2 = np.array(forw) # time x reaction x particle
back2 = np.array(back)
conv2 = np.array(conv)

print(conv2[20,0:50,1])

#for i_batch in range(20):
#    for i_part in [1]:    
#        (figure, axes) = mpl_helper.make_fig(right_margin=1.2)
#        for i_spec in [i_batch*6+0, i_batch*6+1, i_batch*6+2, i_batch*6+3, i_batch*6+4, i_batch*6+5]:
#            print(i_batch,i_part,i_spec)
#            axes.semilogy(time, forw[:,i_spec,i_part]-back[:,i_spec,i_part], label=i_spec)
#            axes.set_xlabel(r"time / min")
#            axes.set_ylabel(r"$R_f - R_b$ reaction rate")
#            axes.grid(True)
#            axes.legend(loc=(1.05,0))
#            figure.savefig(f"figs/urban_plume_react_{i_batch}_{i_part}.pdf")

for i_part in [1]:    
    (figure, axes) = mpl_helper.make_fig(left_margin=1,right_margin=1.8)
    for i_spec in [78]:
        print(i_part,i_spec)
        axes.plot(time[50:601], forw1[50:601,i_spec,i_part], label=str(i_spec+1)+" forw base")
        axes.plot(time[50:601], back1[50:601,i_spec,i_part], linestyle='--',label=str(i_spec+1)+" back base")
        axes.plot(time[50:601], forw2[50:601,i_spec,i_part], label=str(i_spec+1)+" forw nh3")
        axes.plot(time[50:601], back2[50:601,i_spec,i_part], linestyle='--',label=str(i_spec+1)+" back nh3")
        axes.set_xlabel(r"time / min")
        axes.set_ylabel(r"reaction rate")
        axes.grid(True)
        axes.legend(loc=(1.05,0))
        figure.savefig(f"figs/urban_plume_react_{i_part}_{i_spec}.pdf")

for i_part in [1]:    
    (figure, axes) = mpl_helper.make_fig(left_margin=1,right_margin=1.8)
    for i_spec in [146]:
        print(i_part,i_spec)
        axes.plot(time[50:601], forw1[50:601,i_spec,i_part], label=str(i_spec+1)+" forw base")
        axes.plot(time[50:601], back1[50:601,i_spec,i_part], linestyle='--',label=str(i_spec+1)+" back base")
        axes.plot(time[50:601], forw2[50:601,i_spec,i_part], label=str(i_spec+1)+" forw nh3")
        axes.plot(time[50:601], back2[50:601,i_spec,i_part], linestyle='--',label=str(i_spec+1)+" back nh3")
        axes.set_xlabel(r"time / min")
        axes.set_ylabel(r"reaction rate")
        axes.grid(True)
        axes.legend(loc=(1.05,0))
        figure.savefig(f"figs/urban_plume_react_{i_part}_{i_spec}.pdf")

for i_part in [1]:    
    (figure, axes) = mpl_helper.make_fig(left_margin=1,right_margin=1.8)
    for i_spec in [145]:
        print(i_part,i_spec)
        axes.plot(time[50:601], forw1[50:601,i_spec,i_part], label=str(i_spec+1)+" forw base")
        axes.plot(time[50:601], back1[50:601,i_spec,i_part], linestyle='--',label=str(i_spec+1)+" back base")
        axes.plot(time[50:601], forw2[50:601,i_spec,i_part], label=str(i_spec+1)+" forw nh3")
        axes.plot(time[50:601], back2[50:601,i_spec,i_part], linestyle='--',label=str(i_spec+1)+" back nh3")
        axes.set_xlabel(r"time / min")
        axes.set_ylabel(r"reaction rate")
        axes.grid(True)
        axes.legend(loc=(1.05,0))
        figure.savefig(f"figs/urban_plume_react_{i_part}_{i_spec}.pdf")

for i_part in [1]:    
    (figure, axes) = mpl_helper.make_fig(left_margin=1,right_margin=1.8)
    for i_spec in [147]:
        print(i_part,i_spec)
        axes.plot(time[50:601], forw1[50:601,i_spec,i_part], label=str(i_spec+1)+" forw base")
        axes.plot(time[50:601], back1[50:601,i_spec,i_part], linestyle='--',label=str(i_spec+1)+" back base")
        axes.plot(time[50:601], forw2[50:601,i_spec,i_part], label=str(i_spec+1)+" forw nh3")
        axes.plot(time[50:601], back2[50:601,i_spec,i_part], linestyle='--',label=str(i_spec+1)+" back nh3")
        axes.set_xlabel(r"time / min")
        axes.set_ylabel(r"reaction rate")
        axes.grid(True)
        axes.legend(loc=(1.05,0))
        figure.savefig(f"figs/urban_plume_react_{i_part}_{i_spec}.pdf")        
