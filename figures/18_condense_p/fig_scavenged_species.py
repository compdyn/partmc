#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import mpl_helper
import partmc



out_filename1 = "figs/scavenged_bc.pdf" 
out_filename2 = "figs/scavenged_oc.pdf" 
out_filename3 = "figs/scavenged_so4.pdf" 
out_filename4 = "figs/scavenged_no3.pdf" 
out_filename5 = "figs/scavenged_nh4.pdf" 
out_filename6 = "figs/scavenged_num.pdf" 

time_array = np.linspace(0,3,4)
frac_array = np.zeros((2,6,len(time_array))) # run x species x time
frac_num_array = np.zeros((2,5,len(time_array))) # run x species x time
run_list = ["_ne", "_we"]

print len(time_array)

for k in range(0,2):
    run = run_list[k]
    print 'case ', run
    for counter in range(0,4):
        print 'counter ', counter
        if (counter == 0):
            plume_index = 4
        elif (counter == 1):
            plume_index = 8
        elif (counter ==2):
            plume_index = 13
        elif (counter == 3):
            plume_index = 20

        in_dir = "../../scenarios/8_condense_p/out%s/out_%02d/" % (run, plume_index)
        in_filename = "condense_0001_00000061.nc"
        print in_dir+in_filename

        ncf = scipy.io.netcdf.netcdf_file(in_dir+in_filename, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        ncf.close()

        final_wet_diameters = particles.diameters()
        is_activated = (final_wet_diameters > 2e-6)
        total_number = len(particles.ids)
        number_act = sum(is_activated)

        bc = particles.masses(include = ["BC"]) 
        oc = particles.masses(include = ["OC"])
        so4 = particles.masses(include = ["SO4"])
        no3 = particles.masses(include = ["NO3"])
        nh4 = particles.masses(include = ["NH4"])

        contains_bc = ((bc > 0) & (is_activated))
        contains_oc = ((oc > 0) & (is_activated))
        contains_so4 = ((so4 > 0) & (is_activated))
        contains_no3 = ((no3 > 0) & (is_activated))
        contains_nh4 = ((nh4 > 0) & (is_activated))

        print "particles that are activated and contain bc ", sum(contains_bc)

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

        bc_num_act = sum(contains_bc)
        oc_num_act = sum(contains_oc)
        so4_num_act = sum(contains_so4)
        no3_num_act = sum(contains_no3)
        nh4_num_act = sum(contains_nh4)

        frac_array[k,0,counter] = bc_mass_act / total_bc
        frac_array[k,1,counter] = oc_mass_act / total_oc
        frac_array[k,2,counter] = so4_mass_act / total_so4
        frac_array[k,3,counter] = no3_mass_act / total_no3
        frac_array[k,4,counter] = nh4_mass_act / total_nh4
        frac_array[k,5,counter] = float(number_act) / total_number

        frac_num_array[k,0,counter] = float(bc_num_act) / total_number
        frac_num_array[k,1,counter] = float(oc_num_act) / total_number
        frac_num_array[k,2,counter] = float(so4_num_act) / total_number
        frac_num_array[k,3,counter] = float(no3_num_act) / total_number
        frac_num_array[k,4,counter] = float(nh4_num_act) / total_number

        print in_filename, bc_mass_act, total_bc, number_act, total_number 
        print bc_num_act, oc_num_act, so4_num_act, no3_num_act, nh4_num_act
        print bc_num_act, total_number, bc_num_act / total_number

#plt.clf()
#plt.plot(time_array, frac_array[0,5,:], 'kx-', label = 'num ne')
#plt.plot(time_array, frac_array[1,5,:], 'rx-', label = 'num we')
#plt.axis([0,3,0,1])
#plt.legend(loc = 'lower right')
#plt.xlabel("time ")
#plt.ylabel("scavenged fraction")
#plt.grid(True)
#fig = plt.gcf()
#fig.savefig(out_filename6)

plt.clf()
plt.plot(time_array+1, frac_num_array[0,0,:], 'kx', label = 'bc ne')
plt.plot(time_array+1, frac_num_array[1,0,:], 'rx', label = 'bc we')
plt.axis([0,4,0,1])
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
plt.grid(True)
fig = plt.gcf()
fig.savefig(out_filename1)

plt.clf()
plt.plot(time_array, frac_num_array[0,1,:], 'kx', label = 'oc ne')
plt.plot(time_array, frac_num_array[1,1,:], 'rx', label = 'oc we')
plt.axis([0,3,0,1])
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
plt.axis([0,3,0,1])
plt.grid(True)
fig = plt.gcf()
fig.savefig(out_filename2)

plt.clf()
plt.plot(time_array, frac_num_array[0,2,:], 'kx', label = 'so4 ne')
plt.plot(time_array, frac_num_array[1,2,:], 'rx', label = 'so4 we')
plt.axis([0,3,0,1])
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
plt.axis([0,3,0,1])
plt.grid(True)
fig = plt.gcf()
fig.savefig(out_filename3)

plt.clf()
plt.plot(time_array, frac_num_array[0,3,:], 'k-', label = 'no3 ne')
plt.plot(time_array, frac_num_array[1,3,:], 'kx', label = 'no3 we')
plt.axis([0,3,0,1])
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
plt.axis([0,3,0,1])
plt.grid(True)
fig = plt.gcf()
fig.savefig(out_filename4)

plt.clf()
plt.plot(time_array, frac_num_array[0,4,:], 'kx', label = 'nh4 ne')
plt.plot(time_array, frac_num_array[1,4,:], 'mx', label = 'nh4 we')
plt.axis([0,3,0,1])
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
plt.axis([0,3,0,1])
plt.grid(True)
fig = plt.gcf()
fig.savefig(out_filename5)
