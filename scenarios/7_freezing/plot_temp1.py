#!/usr/bin/env python
import os
import glob
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace

sim_days = 2
log_show = False
#draw_ice_type = "num_conc"
#OutDir = "out_base"
#OutDir = "out_coag"
#OutDir = "out_cold"
#OutDir = "out_poster"
#OutDir = "out_poster2"
OutDir = "."
fontsize_label = 15
caseName = "out_case2"

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

t_list = np.array(list(range(1, 864 *sim_days)))
isFirst = True 
for t in t_list:
    fileName = OutDir + "/" + caseName+ "/freezing_part_0001_" + str(t).zfill(8) + ".nc"
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

fig = plt.figure(figsize = (12, 6))
ax = fig.add_subplot(1, 1, 1)
ax.plot(t_list * 100 / 3600, mixing_ratio_ice * 1000, label = "Ice crystals", color = "blue")
ax.set_ylabel("Mixing ratio of ice (g/Kg)", color = "blue", fontsize = 15)
axt = ax.twinx()
axt.plot(t_list * 100 / 3600, temperature_data - 273.15, label = "temperature", color = "red")
axt.set_ylabel("Temperature (˚C)", color = "red", fontsize = 15)
axt.plot(t_list * 100 / 3600, np.full((len(t_list)), 0), linestyle = "dashed", color = "grey")
ax.set_xlabel("time (hour)", fontsize = 15)
ax.grid()
#ax.legend()
plt.show()


