#!/usr/bin/env python2.5

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

in_filename = "../../scenarios/1_urban_plume/out/urban_plume_wc_0001_00000025.nc"
out_filename = "figs/ccn_spectrum_wc_24.pdf"
ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
particles = partmc.aero_particle_array_t(ncf)
env_state = partmc.env_state_t(ncf)
ncf.close()

s_crit = (particles.critical_rel_humids(env_state) - 1)*100
total_number = len(particles.ids)
ss_array = np.logspace(-2,1,100)
frac_array = np.zeros_like(ss_array)

for i in range(len(ss_array)):
    activated = (s_crit < ss_array[i])
    number_act = sum(activated)
    frac_array[i] = float(number_act) / total_number

plt.grid(True)
plt.semilogx(ss_array, frac_array)
plt.ylabel("CCN/CN")
plt.xlabel("critical supersaturation (%)")
fig = plt.gcf()
fig.savefig(out_filename)

