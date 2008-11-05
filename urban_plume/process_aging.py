#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
from fig_helper import *

time_filename_list_nc = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
env_state = read_any(env_state_t, netcdf_dir_nc, netcdf_pattern_nc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, filename, key] in time_filename_list_nc]) / 60

outf_a = open("aging_a_nc", "w")
outf_p = open("aging_p_nc", "w")
outf_ea = open("aging_ea_nc", "w")
outf_ep = open("aging_ep_nc", "w")

old_id = set()

for [time, filename, key] in time_filename_list_nc:
    print time, filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()
    num_den = 1.0 / particles.comp_vol
    critical_ss = particles.kappa_rh(env_state) - 1.0
    total_num_den = num_den.sum()
    outf_a.write("%f " % time)
    outf_p.write("%f " % time)
    for ss_activ in [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01]:
        aged_num_den = 0.0
        fresh_num_den = 0.0
        aged_emissions = 0.0
        fresh_emissions = 0.0
        
        for i in range(num_den.size):
            # i is particle index
            if particles.id[i] not in old_id:    
                if critical_ss[i] < ss_activ:
                    aged_emissions += num_den[i]
                else:
                    fresh_emissions += num_den[i]
            if critical_ss[i] < ss_activ:
                aged_num_den += num_den[i]
            else:
                fresh_num_den += num_den[i]
        print time / 60., ss_activ, aged_emissions, fresh_emissions, aged_num_den, fresh_num_den, total_num_den
        outf_a.write("%e " % aged_num_den)
        outf_p.write("%e " % fresh_num_den)
        outf_ea.write("%e " % aged_emissions)
        outf_ep.write("%e " % fresh_emissions)
    outf_a.write("\n")
    outf_p.write("\n")
    outf_ea.write("\n")
    outf_ep.write("\n")
    old_id = set()
    for i in range(num_den.size):
        old_id.add(particles.id[i]) 
outf_a.close()
outf_p.close()
outf_ea.close()
outf_ep.close()

