#!/usr/bin/env python
# Copyright (C) 2007-2009 Nicole Riemer and Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
sys.path.append("../tool")
from pmc_data_nc import *
from fig_helper import *

time_filename_list = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
env_state = read_any(env_state_t, netcdf_dir_nc, netcdf_pattern_nc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, filename, key] in time_filename_list]) / 60

outf_a = open("out/aging_aged_nc", "w")
outf_f = open("out/aging_fresh_nc", "w")
outf_ea = open("out/aging_emission_aged_nc", "w")
outf_ef = open("out/aging_emission_fresh_nc", "w")
outf_la = open("out/aging_loss_aged_nc", "w")
outf_lf = open("out/aging_loss_fresh_nc", "w")
outf_ta = open("out/aging_transfer_to_aged_nc", "w")
outf_tf = open("out/aging_transfer_to_fresh_nc", "w")
outf_h = open("out/aging_height_nc", "w")

old_id = set()

const = load_constants("../src/constants.f90")

first_time = True
for [time, filename, key] in time_filename_list:
    print time, filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    particles.id = [int(i) for i in particles.id]
    env_state = env_state_t(ncf)
    ncf.close()
    num_den = 1.0 / array(particles.comp_vol)
    soot_mass = particles.mass(include = ["BC"])
    critical_ss = particles.kappa_rh(env_state, const) - 1.0
    total_num_den = num_den.sum()
    outf_h.write("%f %e\n" % (time, env_state.height))
    outf_a.write("%f " % time)
    outf_f.write("%f " % time)
    outf_ea.write("%f " % time)
    outf_ef.write("%f " % time)
    outf_la.write("%f " % time)
    outf_lf.write("%f " % time)
    outf_ta.write("%f " % time)
    outf_tf.write("%f " % time)

    for ss_activ in [0.001, 0.003, 0.006, 0.010]:
        aged_num_den = 0.0
        fresh_num_den = 0.0
        aged_emissions = 0.0
        fresh_emissions = 0.0
        aged_loss = 0.0
        fresh_loss = 0.0
        transfer_to_aged = 0.0
        transfer_to_fresh = 0.0

        current_id = set()
        for i in range(particles.n_particles):
            current_id.add(particles.id[i])

        if not first_time:
            old_id = set()
            for i in range(old_particles.n_particles):
                old_id.add(old_particles.id[i])

            old_aged_id = set()
            old_fresh_id = set()
            for i in range(old_particles.n_particles):
                if old_critical_ss[i] < ss_activ:
                    old_aged_id.add(old_particles.id[i])
                else:
                    old_fresh_id.add(old_particles.id[i])

            for i in range(old_particles.n_particles):
                if old_soot_mass[i] > 0.0:
                    if old_particles.id[i] not in current_id:
                        if old_critical_ss[i] < ss_activ:
                            aged_loss += old_num_den[i]
                        else:
                            fresh_loss += old_num_den[i]

        for i in range(particles.n_particles):
            if soot_mass[i] > 0.0:
                if not first_time:
                    if particles.id[i] not in old_id:    
                        if critical_ss[i] < ss_activ:
                            aged_emissions += num_den[i]
                        else:
                            fresh_emissions += num_den[i]
                if critical_ss[i] < ss_activ:
                    aged_num_den += num_den[i]
                    if particles.id[i] in old_fresh_id:
                        transfer_to_aged += num_den[i]
                else:
                    fresh_num_den += num_den[i]
                    if particles.id[i] in old_aged_id:
                        transfer_to_fresh += num_den[i]
        print "%e %e %e %e %e %e %e %e %e %e %e" \
            % (time / 60., ss_activ, aged_emissions,
               fresh_emissions, aged_loss, fresh_loss,
               aged_num_den, fresh_num_den, total_num_den,
               transfer_to_aged, transfer_to_fresh)
        outf_a.write("%e " % aged_num_den)
        outf_f.write("%e " % fresh_num_den)
        outf_ea.write("%e " % aged_emissions)
        outf_ef.write("%e " % fresh_emissions)
        outf_la.write("%e " % aged_loss)
        outf_lf.write("%e " % fresh_loss)
        outf_ta.write("%e " % transfer_to_aged)
        outf_tf.write("%e " % transfer_to_fresh)
    outf_a.write("\n")
    outf_f.write("\n")
    outf_ea.write("\n")
    outf_ef.write("\n")
    outf_la.write("\n")
    outf_lf.write("\n")
    outf_ta.write("\n")
    outf_tf.write("\n")
    old_particles = particles
    old_num_den = num_den
    old_critical_ss = critical_ss
    old_soot_mass = soot_mass
    first_time = False

outf_a.close()
outf_f.close()
outf_ea.close()
outf_ef.close()
outf_la.close()
outf_lf.close()
outf_ta.close()
outf_tf.close()
outf_h.close()
