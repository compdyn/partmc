#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import scipy
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
const = pmc_data_nc.load_constants("../../src/constants.f90")

def effective_diam(env_state, const, kappa, s_crit):
    def f(dry_diam):
        return pmc_data_nc.critical_rh(env_state, const, np.array([kappa]), np.array([dry_diam]))[0] - (s_crit + 1)
    return scipy.optimize.brentq(f, 1e-10, 1)

in_filename = "../../urban_plume/out/urban_plume_wc_0001_00000007.nc"
out_filename = "figs/ccn_spectrum_diam.pdf"
ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename)
particles = pmc_data_nc.aero_particle_array_t(ncf)
env_state = pmc_data_nc.env_state_t(ncf)
ncf.close()

s_crit = (particles.critical_rh(env_state, const) - 1)*100
kappa = particles.solute_kappa()
dry_volumes = particles.volume(exclude = ["H2O"]) 
average_kappa = 1/sum(dry_volumes) * sum(kappa * dry_volumes)
print average_kappa
dry_diameter = particles.dry_diameter()
total_number = len(particles.id)
ss_array = np.logspace(-2,1,100)
frac_array_true = np.zeros_like(ss_array)
frac_array_diam = np.zeros_like(ss_array)
for i in range(len(ss_array)):
    d_eff = effective_diam(env_state, const, average_kappa, ss_array[i]/100)
    print ss_array[i], d_eff
    activated_true = (s_crit < ss_array[i])
    number_act_true = sum(activated_true)
    frac_array_true[i] = float(number_act_true) / total_number

    activated_diam = (dry_diameter > d_eff)
    number_act_diam = sum(activated_diam)
    frac_array_diam[i] = float(number_act_diam)/ total_number

plt.grid(True)
plt.semilogx(ss_array, frac_array_true, label = "true")
plt.semilogx(ss_array, frac_array_diam, label = "diam")
plt.legend()
plt.ylabel("CCN/CN")
plt.xlabel("critical supersaturation (%)")
fig = plt.gcf()
fig.savefig(out_filename)

