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

in_dir = "../../scenarios/3_condense/out/"
out_filename1 = "figs/scavenged_bcwc.pdf" 
out_filename2 = "figs/scavenged_oc_wc.pdf" 
out_filename3 = "figs/scavenged_so4_wc.pdf" 
out_filename4 = "figs/scavenged_no3_wc.pdf" 
out_filename5 = "figs/scavenged_nh4_wc.pdf" 

time_array = np.linspace(0,48,49)
frac_array = np.zeros((4,5,len(time_array))) # run x species x time
run_list = ["ref", "comp", "size", "both"]

for k in range(0,4):
    run = run_list[k]
    for counter in range(1,41):
        in_filename = "cond_%02d_%s_0001_00000601.nc" % (counter, run)
        ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename)
        particles = partmc.aero_particle_array_t(ncf)
        ncf.close()

        final_wet_diameters = particles.diameters()
        is_activated = (final_wet_diameters > 3e-6)

        bc = particles.masses(include = ["BC"]) 
        oc = particles.masses(include = ["OC"])
        so4 = particles.masses(include = ["SO4"])
        no3 = particles.masses(include = ["NO3"])
        nh4 = particles.masses(include = ["NH4"])

        total_bc = sum(bc)
        total_oc = sum(oc)
        total_so4 = sum(so4)
        total_no3 = sum(no3)
        total_nh4 = sum(nh4)

        bc_mass_act = sum(bc[is_activated])
        oc_mass_act = sum(oc[is_activated])
        so4_mass_act = sum(so4[is_activated])
        no3_mass_act = sum(no3[is_activated])
        nh4_mass_act = sum(nh4[is_activated])

        frac_array[k,0,counter] = bc_mass_act / total_bc
        frac_array[k,1,counter] = oc_mass_act / total_oc
        frac_array[k,2,counter] = so4_mass_act / total_so4
        frac_array[k,3,counter] = no3_mass_act / total_no3
        frac_array[k,4,counter] = nh4_mass_act / total_nh4

        print 'bc ', run, counter, total_bc, bc_mass_act

plt.clf()
plt.plot(time_array, frac_array[0,0,:], 'r-', label = 'bc ref')
plt.plot(time_array, frac_array[1,0,:], 'r.', label = 'bc comp')
plt.plot(time_array, frac_array[2,0,:], 'r--', label = 'bc size')
plt.plot(time_array, frac_array[3,0,:], 'r-.', label = 'bc both')
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
fig = plt.gcf()
fig.savefig(out_filename1)

plt.clf()
plt.plot(time_array, frac_array[0,1,:], 'b-', label = 'oc ref')
plt.plot(time_array, frac_array[1,1,:], 'b.', label = 'oc comp')
plt.plot(time_array, frac_array[2,1,:], 'b--', label = 'oc size')
plt.plot(time_array, frac_array[3,1,:], 'b-.', label = 'oc both')
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
fig = plt.gcf()
fig.savefig(out_filename2)

plt.clf()
plt.plot(time_array, frac_array[0,2,:], 'g-', label = 'so4 ref')
plt.plot(time_array, frac_array[1,2,:], 'g.', label = 'so4 comp')
plt.plot(time_array, frac_array[2,2,:], 'g--', label = 'so4 both')
plt.plot(time_array, frac_array[3,2,:], 'g-.', label = 'so4 ref')
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
fig = plt.gcf()
fig.savefig(out_filename3)

plt.clf()
plt.plot(time_array, frac_array[0,3,:], 'k-', label = 'no3 ref')
plt.plot(time_array, frac_array[1,3,:], 'k.', label = 'no3 comp')
plt.plot(time_array, frac_array[2,3,:], 'k--', label = 'no3 size')
plt.plot(time_array, frac_array[3,3,:], 'k-.', label = 'no3 both')
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
fig = plt.gcf()
fig.savefig(out_filename4)

plt.clf()
plt.plot(time_array, frac_array[0,4,:], 'm-', label = 'nh4 ref')
plt.plot(time_array, frac_array[1,4,:], 'm.', label = 'nh4 comp')
plt.plot(time_array, frac_array[2,4,:], 'm--', label = 'nh4 size')
plt.plot(time_array, frac_array[3,4,:], 'm-.', label = 'nh4 both')
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
fig = plt.gcf()
fig.savefig(out_filename5)
