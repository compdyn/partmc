#!/usr/bin/env python2.5

import scipy.io
import sys
import scipy
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def effective_diam(env_state, kappa, s_crit):
    def f(dry_diam):
        return partmc.critical_rel_humids(env_state, np.array([kappa]), np.array([dry_diam]))[0] - (s_crit + 1)
    return scipy.optimize.brentq(f, 1e-10, 1)

in_filename = "../../scenarios/1_urban_plume/out/urban_plume_wc_0001_00000013.nc"
out_filename = "figs/ccn_spectrum_diam.pdf"
ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
particles = partmc.aero_particle_array_t(ncf)
env_state = partmc.env_state_t(ncf)
ncf.close()

s_crit = (particles.critical_rel_humids(env_state) - 1)*100
bc = particles.masses(include = ["BC"])
total_bc = sum(bc)

kappa = particles.kappas()
dry_volumes = particles.volumes(exclude = ["H2O"]) 
average_kappa = 1/sum(dry_volumes) * sum(kappa * dry_volumes)
print average_kappa
dry_diameters = particles.dry_diameters()
total_number = len(particles.ids)
ss_array = np.logspace(-2,1,100)
frac_array_true = np.zeros_like(ss_array)
frac_array_diam = np.zeros_like(ss_array)
for i in range(len(ss_array)):
    d_eff = effective_diam(env_state, average_kappa, ss_array[i]/100)
    print ss_array[i], d_eff
    activated_true = (s_crit < ss_array[i])
    bc_mass_act_true = sum(bc[activated_true])
    frac_array_true[i] = bc_mass_act_true / total_bc

    activated_diam = (dry_diameters > d_eff)
    bc_mass_act_diam = sum(bc[activated_diam])
    frac_array_diam[i] = bc_mass_act_diam/ total_bc

plt.grid(True)
plt.semilogx(ss_array, frac_array_true, label = "true")
plt.semilogx(ss_array, frac_array_diam, label = "diam")
plt.legend()
plt.ylabel("scavenged BC mass")
plt.xlabel("critical supersaturation (%)")
fig = plt.gcf()
fig.savefig(out_filename)

