#!/usr/bin/env python
import os
import glob
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace


log_show = False
draw_relhum = False
#draw_ice_type = "num_conc"
draw_ice_type = "mix"
#OutDir = "out_base"
#OutDir = "out_coag"
#OutDir = "out_cold"
#OutDir = "out_poster"
#OutDir = "out_case1"
OutDir = "output/evaporation_base"

fontsize_label = 15
assert draw_ice_type in ["num_conc", "mix"]
t_max = 26
t_output = 1
nFile = int(t_max / t_output) + 1

temperature_data = []
pressure_data = []
relhum_data = []
freezing_data = []
P_frozen_data = []
num_conc_data = []
freezing_total = []
freezing_total2 = []
freezing_mass_data = [] 
freezing_mass_data2 = [] 

#t_list = list(range(1, 864 *2))
id_file_list = np.array(list(range(1, nFile + 1)))
t_list = (id_file_list - 1) * t_output
isFirst = True
for id_file in id_file_list:
    fileName = OutDir + "/freezing_part_0001_" + str(id_file).zfill(8) + ".nc"
    ncf = nc.Dataset(fileName)
    part_frozen = ncf.variables["aero_frozen"][:].filled(np.nan)
    p_frozen = ncf.variables["aero_frozen_probability"][:].filled(np.nan)
    num_conc = ncf.variables["aero_num_conc"][:].filled(np.nan)
    particle_mass = ncf.variables["aero_particle_mass"][:].filled(np.nan).sum(axis = 0)
    temperature = float(ncf.variables["temperature"][0].filled(np.nan))
    relhum = float(ncf.variables["relative_humidity"][0].filled(np.nan))
    pressure = float(ncf.variables["pressure"][0].filled(np.nan))

    freezing_data.append(part_frozen)
    P_frozen_data.append(p_frozen)
    num_conc_data.append(num_conc)
    freezing_total.append(int(np.sum(p_frozen * num_conc)))
    freezing_total2.append(int(np.sum(part_frozen * num_conc)))
    freezing_mass_data.append( np.sum(p_frozen * num_conc * particle_mass) )
    freezing_mass_data2.append( np.sum(part_frozen * num_conc * particle_mass) )

    temperature_data.append(temperature)
    relhum_data.append(relhum * 100)
    pressure_data.append(pressure)


"""
freezing_data = np.array(freezing_data)
P_frozen_data = np.array(P_frozen_data)
num_conc_data = np.array(num_conc_data)
freezing_total = np.array(freezing_total)
freezing_total2 = np.array(freezing_total2)
"""
temperature_data = np.array(temperature_data)
relhum_data = np.array(relhum_data)
pressure_data = np.array(pressure_data)
freezing_mass_data = np.array(freezing_mass_data)
freezing_mass_data2 = np.array(freezing_mass_data2)

air_density_data = pressure / (286 * temperature_data)
mixing_ratio_ice = freezing_mass_data / (air_density_data * 1)
mixing_ratio_ice2 = freezing_mass_data2 / (air_density_data * 1)

f_ufz_mean = np.array([np.mean(1 - p_frozen) for p_frozen in P_frozen_data])
f_ufz_std = np.array([np.std(1 - p_frozen) for p_frozen in P_frozen_data])
P_frz_mean = np.array([np.mean(p_frozen) for p_frozen in P_frozen_data])
P_frz_std = np.array([np.std(p_frozen) for p_frozen in P_frozen_data])

#freezing_data = np.array(freezing_data)
#set_trace()
#freezing_total = freezing_data.sum(axis = 1)

#fig = plt.figure(figsize = (12, 10))
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (12, 10), sharex = True)
#ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(t_list, temperature_data, color = "red")
ax1.plot(t_list, np.full((len(t_list)), 273.15), color = "red" if draw_relhum else "blue", linestyle = "dashed")
#ax1.set_xlabel("time (s)")
ax1.set_ylabel("Temperature (K)", color = "red" if draw_relhum else "black", fontsize = fontsize_label)
ax1.grid()
if draw_relhum:
    ax1t = ax1.twinx()
    ax1t.plot(t_list, relhum_data, color = "blue")
    ax1t.set_ylabel("Relative humidity (%)", color = "blue", fontsize = fontsize_label)

#ax2=fig.add_subplot(3, 1, 2)
if draw_ice_type == "num_conc":
    #ax2.plot(t_list, np.log10(freezing_total) if log_show else freezing_total, label = "probability", color = "orange")
    ax2.plot(t_list, np.log10(freezing_total2) if log_show else freezing_total2, label = "true false", color = "blue")
    ax2.set_ylabel("frozen concentration " + ("10^* (m^-3)" if log_show else "(m^-3)"), fontsize = fontsize_label)
else:
    #ax2.plot(t_list, mixing_ratio_ice * 1000, label = "probability", color = "orange")
    ax2.plot(t_list, mixing_ratio_ice2 * 1000, label = "true false", color = "blue")
    ax2.set_ylabel("ice mixing ratio (g/kg)", fontsize = fontsize_label)
ax2.grid()
#ax2.set_xlim([t_list[0], t_list[-1]])
#ax2.set_xlabel("time (s)")

ax2.legend()

#ax3 = fig.add_subplot(3, 1, 3)
if log_show:
    ax3.plot(t_list, np.log10(f_ufz_mean), color = "red")
    ax3.fill_between(t_list, np.log10(f_ufz_mean - f_ufz_std), np.log10(f_ufz_mean + f_ufz_std), alpha = 0.1, color = "red", edgecolor = None)

else:
    #ax3.plot(t_list, f_ufz_mean, color = "red")
    #ax3.fill_between(t_list, f_ufz_mean - f_ufz_std, f_ufz_mean + f_ufz_std, alpha = 0.1, color = "red", edgecolor = None)
    ax3.plot(t_list, P_frz_mean, color = "red")
    ax3.fill_between(t_list, P_frz_mean - P_frz_std, P_frz_mean + P_frz_std, alpha = 0.1, color = "red", edgecolor = None)

    ax3.set_ylim([-0.05, 1.05])
#ax2.set_xlim([t_list[0], t_list[-1]])
ax3.grid()
ax3.set_xlabel("time (s)", fontsize = fontsize_label)
ax3.set_ylabel("frozen ratio" + (" (10^*)" if log_show else ""), fontsize = fontsize_label)

fig.subplots_adjust(hspace = 0.1)
fig.suptitle("Immersion freezing (partMC-ABIFM)", fontsize = 20)

plt.show()

