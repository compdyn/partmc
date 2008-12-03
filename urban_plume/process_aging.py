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

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, filename, key] in time_filename_list]) / 60

outf_a = open("out/aging_a_wc", "w")
outf_p = open("out/aging_p_wc", "w")
outf_ea = open("out/aging_ea_wc", "w")
outf_ep = open("out/aging_ep_wc", "w")
outf_h = open("out/aging_h_wc", "w")

old_id = set()

const = load_constants("../src/constants.f90")

for [time, filename, key] in time_filename_list:
    print time, filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()
    num_den = 1.0 / particles.comp_vol
    critical_ss = particles.kappa_rh(env_state, const) - 1.0
    total_num_den = num_den.sum()
    outf_h.write("%f %e\n" % (time, env_state.height))
    outf_a.write("%f " % time)
    outf_f.write("%f " % time)
    outf_ea.write("%f " % time)
    outf_ef.write("%f " % time)

    for ss_activ in [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01]:
        aged_num_den = 0.0
        fresh_num_den = 0.0
        aged_emissions = 0.0
        fresh_emissions = 0.0
        
        for i in range(particles.n_particles):
            # i is particle index
            print particles.id[i]
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
        outf_f.write("%e " % fresh_num_den)
        outf_ea.write("%e " % aged_emissions)
        outf_ef.write("%e " % fresh_emissions)
    outf_a.write("\n")
    outf_f.write("\n")
    outf_ea.write("\n")
    outf_ef.write("\n")
    old_id = set()
    for i in range(num_den.size):
        old_id.add(particles.id[i]) 
outf_a.close()
outf_f.close()
outf_ea.close()
outf_ef.close()
outf_h.close()

