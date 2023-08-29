#!/usr/bin/env python
import os
import glob
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace

#sim_days = 2
log_show = False
#draw_ice_type = "num_conc"
draw_ice_type = "mix"
#OutDir = "out_base"
#OutDir = "out_coag"
#OutDir = "out_cold"
#OutDir = "out_poster"
#OutDir = "out_poster2"
OutDir = "./output"
fontsize_label = 15
assert draw_ice_type in ["num_conc", "mix"]
#casesName = ["out_case1", "out_case2", "out_case3"]
#casesName = ["out_case1", "out_case2", "out_case3"]
#casesName = ["out_new2", "out_new3", "out_new4", "out_old2", "out_old3", "out_old4"]
#casesLabel = ["new", "new", "new", "old", "old", "old"]
#casesColor = ["red", "red", "red",  "blue", "blue", "blue"]

casesName = ["mix_OIN100", "mix_ILT100", "mix_OIN50_ILT50"]
casesLabel = ["mix_OIN100", "mix_ILT100", "mix_OIN50_ILT50"]
#casesLabel = ["condense_t1", "test_m2"]
casesColor = ["red", "blue", "green"]

t_max = 3600 * 2
t_output = 100
nFile = int(t_max / t_output) + 1


def plot_case(caseName, caseLabel, caseColor, ax1, ax2, ax3):
    print("Processing " + caseName + " ...")
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

    #t_list = list(range(1, 864 *sim_days))
    id_file_list = np.array(list(range(1, nFile + 1)))
    t_list = (id_file_list - 1) * t_output

    isFirst = True
    #for t in t_list:
    for id_file in id_file_list:
        fileName = OutDir + "/" + caseName+ "/freezing_part_0001_" + str(id_file).zfill(8) + ".nc"
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
    
    ax1.plot(t_list, temperature_data, color = caseColor, label = caseLabel)
    ax1.plot(t_list, np.full((len(t_list)), 273.15), color = "grey" , linestyle = "dashed")

    if draw_ice_type == "num_conc":
        #ax2.plot(t_list, np.log10(freezing_total) if log_show else freezing_total, label = caseLabel, color = caseColor)
        ax2.plot(t_list, np.log10(freezing_total2) if log_show else freezing_total2, label = caseLabel, color = caseColor)

    else:
        #ax2.plot(t_list, mixing_ratio_ice * 1000, label = caseLabel, color = caseColor)
        ax2.plot(t_list, mixing_ratio_ice2 * 1000, label = caseLabel, color = caseColor)
    if log_show:
        ax3.plot(t_list, np.log10(f_ufz_mean), color = caseColor, label = caseLabel)
        ax3.fill_between(t_list, np.log10(f_ufz_mean - f_ufz_std), np.log10(f_ufz_mean + f_ufz_std), alpha = 0.1, color = caseColor, edgecolor = None)

    else:
        ax3.plot(t_list, P_frz_mean, color = caseColor, label = caseLabel)
        ax3.fill_between(t_list, P_frz_mean - P_frz_std, P_frz_mean + P_frz_std, alpha = 0.1, color = caseColor, edgecolor = None)

if __name__ == "__main__":

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize = (12, 10), sharex = True)

    for caseName, caseLabel, caseColor in zip(casesName, casesLabel, casesColor):
        plot_case(caseName, caseLabel, caseColor, ax1, ax2, ax3)

    ax1.set_ylabel("temperature (K)", color = "black", fontsize = fontsize_label)
    ax1.grid()
    ax1.legend()

    ax2.set_ylabel("frozen concentration " + ("10^* (m^-3)" if log_show else "(m^-3)"), fontsize = fontsize_label)
    ax2.set_ylabel("ice mixing ratio (g/kg)", fontsize = fontsize_label)
    ax2.grid()

    ax2.legend()


    ax3.set_ylim([-0.05, 1.05])
    ax3.grid()
    ax3.set_xlabel("time (s)", fontsize = fontsize_label)
    ax3.set_ylabel("frozen ratio" + (" (10^*)" if log_show else ""), fontsize = fontsize_label)
    ax3.legend()

    fig.subplots_adjust(hspace = 0.1)
    fig.suptitle("Immersion freezing (partMC-ABIFM)", fontsize = 20)

    plt.show()

