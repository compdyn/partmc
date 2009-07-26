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
from numpy import *

const = load_constants("../src/constants.f90")

for coag in [True, False]:
    if coag:
        time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
        env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
        coag_suffix = "wc"
    else:
        time_filename_list = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
        env_state = read_any(env_state_t, netcdf_dir_nc, netcdf_pattern_nc)
        coag_suffix = "nc"
    start_time_of_day_min = env_state.start_time_of_day / 60
    max_time_min = max([time for [time, filename, key] in time_filename_list]) / 60

    first_time = True
    for [time, filename, key] in time_filename_list:
        #DEBUG
        #if time > 121:
        #    sys.exit(0)
        #DEBUG
        print time, filename
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        particles.id = [int(i) for i in particles.id]
        env_state = env_state_t(ncf)
        ncf.close()
        num_den = 1.0 / array(particles.comp_vol)
        total_num_den = num_den.sum()
        soot_mass = particles.mass(include = ["BC"])
        critical_ss = particles.kappa_rh(env_state, const) - 1.0
        ss_bin = ss_active_axis.find_clipped_outer(critical_ss)
        id_set = set([particles.id[i] for i in range(particles.n_particles)])
        particle_index_by_id = dict([[particles.id[i], i] for i in range(particles.n_particles)])

        num = zeros(n_level_bin + 2, int)
        num_emit = zeros(n_level_bin + 2, int)
        num_dilution = zeros(n_level_bin + 2, int)
        num_halving = zeros(n_level_bin + 2, int)
        num_cond = zeros([n_level_bin + 2, n_level_bin + 2], int)
        num_coag_gain = zeros(n_level_bin + 2, int)
        num_coag_loss = zeros([n_level_bin + 2, n_level_bin + 2], int)

        mass = zeros(n_level_bin + 2, float)
        mass_emit = zeros(n_level_bin + 2, float)
        mass_dilution = zeros(n_level_bin + 2, float)
        mass_halving = zeros(n_level_bin + 2, float)
        mass_cond = zeros([n_level_bin + 2, n_level_bin + 2], float)
        mass_coag_gain = zeros(n_level_bin + 2, float)
        mass_coag_loss = zeros([n_level_bin + 2, n_level_bin + 2], float)

        if any(array(particles.aero_removed_action) == AERO_INFO_HALVED):
            halving_occured = True
        else:
            halving_occured = False

        # num, mass
        for i in range(particles.n_particles):
            if soot_mass[i] > 0.0:
                #DEBUG
                #if critical_ss[i] < 0.001:
                #    print critical_ss[i], ss_bin[i], ss_active_axis.find(array([critical_ss[i]])), ss_active_axis.min
                #DEBUG
                num[ss_bin[i]] += 1
                mass[ss_bin[i]] += soot_mass[i]
        #DEBUG
        #sys.exit(0)
        #DEBUG

        if not first_time:
            removed_particles = {}
            for i in range(len(particles.aero_removed_id)):
                removed_particles[particles.aero_removed_id[i]] = [particles.aero_removed_action[i],
                                                                   particles.aero_removed_other_id[i]]

            if not halving_occured:
                if old_id_set - id_set != set(removed_particles.keys()):
                    print "old_id_set - id_set: ", old_id_set - id_set
                    print "removed_particles: ", set(removed_particles.keys())
                    print "xor: ", ((old_id_set - id_set) ^ set(removed_particles.keys()))
                    raise Exception("lost particle mismatch at t = %f" % time)
            #DEBUG
            #print "len(id_set - old_id_set): ", len(id_set - old_id_set)
            #print "id_set - old_id_set: ", id_set - old_id_set
            #DEBUG

            final_outcomes = {}
            for (id, [action, other_id]) in removed_particles.iteritems():
                final_id = id
                while final_id in removed_particles.keys():
                    [final_action, final_other_id] = removed_particles[final_id]
                    if final_action == AERO_INFO_COAG:
                        final_id = final_other_id
                    else:
                        break
                final_outcomes[id] = [final_action, final_other_id]
            for (id, [action, other_id]) in removed_particles.iteritems():
                if other_id not in final_outcomes.keys():
                    final_outcomes[other_id] = [AERO_INFO_COAG, other_id]
            # final_outcomes now stores the final outcome
            # (diluted, halved, coagulated) of all particles that
            # are removed or which underwent
            # coagulation. Particles that are still present but
            # which underwent coagulation will be in
            # final_outcomes with other_id set to id.

            # dilution, halving, coag_loss
            for old_i in range(old_particles.n_particles):
                if old_soot_mass[old_i] > 0.0:
                    old_id = old_particles.id[old_i]
                    if old_id in final_outcomes.keys():
                        [final_action, final_other_id] = final_outcomes[old_id]
                        if final_action == AERO_INFO_NONE:
                            raise Exception("found AERO_INFO_NONE at t = %f" % time)
                        elif final_action == AERO_INFO_DILUTION:
                            num_dilution[old_ss_bin[old_i]] += 1
                            mass_dilution[old_ss_bin[old_i]] += old_soot_mass[old_i]
                        elif final_action == AERO_INFO_COAG:
                            final_i = particle_index_by_id[final_other_id]
                            num_coag_loss[old_ss_bin[old_i], ss_bin[final_i]] += 1
                            mass_coag_loss[old_ss_bin[old_i], ss_bin[final_i]] += old_soot_mass[old_i]
                        elif final_action == AERO_INFO_HALVED:
                            num_halving[old_ss_bin[old_i]] += 1
                            mass_halving[old_ss_bin[old_i]] += old_soot_mass[old_i]

            coag_id = set()
            for (id, [action, other_id]) in removed_particles.iteritems():
                if action == AERO_INFO_COAG:
                    coag_id.add(id)
                    coag_id.add(other_id)

            if not (coag_id <= old_id_set):
                raise Exception("coag_id not a subset of old_id_set at t = %f" % time)
            if not (coag_id <= set(final_outcomes.keys())):
                raise Exception("coag_id not a subset of final_outcomes at t = %f" % time)

            # num_emit, num_cond, num_coag_gain
            for i in range(particles.n_particles):
                id = particles.id[i]
                if soot_mass[i] > 0.0:
                    if id in coag_id:
                        num_coag_gain[ss_bin[i]] += 1
                        mass_coag_gain[ss_bin[i]] += soot_mass[i]
                    elif id in old_id_set:
                        old_i = old_particle_index_by_id[id]
                        num_cond[old_ss_bin[old_i], ss_bin[i]] += 1
                        mass_cond[old_ss_bin[old_i], ss_bin[i]] += soot_mass[i]
                    else:
                        #DEBUG
                        #print "emit id: ", id
                        #DEBUG
                        num_emit[ss_bin[i]] += 1
                        mass_emit[ss_bin[i]] += soot_mass[i]

        time_array = array([time])
        height_array = array([env_state.height])
        comp_vol_array = array([particles.comp_vol[0]])

        filename = os.path.join(aging_data_dir,
                                "aging_%s_%s_%%s.txt" % (coag_suffix, key))

        savetxt(filename % "time", time_array)
        savetxt(filename % "height", height_array)
        savetxt(filename % "comp_vol", comp_vol_array)

        savetxt(filename % "num", num, fmt = "%d")
        savetxt(filename % "num_emit", num_emit, fmt = "%d")
        savetxt(filename % "num_dilution", num_dilution, fmt = "%d")
        savetxt(filename % "num_halving", num_halving, fmt = "%d")
        savetxt(filename % "num_cond", num_cond, fmt = "%d")
        savetxt(filename % "num_coag_gain", num_coag_gain, fmt = "%d")
        savetxt(filename % "num_coag_loss", num_coag_loss, fmt = "%d")

        savetxt(filename % "mass", mass)
        savetxt(filename % "mass_emit", mass_emit)
        savetxt(filename % "mass_dilution", mass_dilution)
        savetxt(filename % "mass_halving", mass_halving)
        savetxt(filename % "mass_cond", mass_cond)
        savetxt(filename % "mass_coag_gain", mass_coag_gain)
        savetxt(filename % "mass_coag_loss", mass_coag_loss)

        old_particles = particles
        old_num_den = num_den
        old_soot_mass = soot_mass
        old_critical_ss = critical_ss
        old_ss_bin = ss_bin
        old_id_set = id_set
        old_particle_index_by_id = particle_index_by_id
        first_time = False
