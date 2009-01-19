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

    outf_height = open("out/aging_%s_height.txt" % coag_suffix, "w")
    outf_comp_vol = open("out/aging_%s_comp_vol.txt" % coag_suffix, "w")

    outf_num_a = open("out/aging_%s_num_a.txt" % coag_suffix, "w")
    outf_num_f = open("out/aging_%s_num_f.txt" % coag_suffix, "w")
    outf_num_emit_a = open("out/aging_%s_num_emit_a.txt" % coag_suffix, "w")
    outf_num_emit_f = open("out/aging_%s_num_emit_f.txt" % coag_suffix, "w")
    outf_num_dilution_a = open("out/aging_%s_num_dilution_a.txt" % coag_suffix, "w")
    outf_num_dilution_f = open("out/aging_%s_num_dilution_f.txt" % coag_suffix, "w")
    outf_num_halving_a = open("out/aging_%s_num_halving_a.txt" % coag_suffix, "w")
    outf_num_halving_f = open("out/aging_%s_num_halving_f.txt" % coag_suffix, "w")
    outf_num_cond_a_a = open("out/aging_%s_num_cond_a_a.txt" % coag_suffix, "w")
    outf_num_cond_a_f = open("out/aging_%s_num_cond_a_f.txt" % coag_suffix, "w")
    outf_num_cond_f_a = open("out/aging_%s_num_cond_f_a.txt" % coag_suffix, "w")
    outf_num_cond_f_f = open("out/aging_%s_num_cond_f_f.txt" % coag_suffix, "w")
    outf_num_coag_gain_a = open("out/aging_%s_num_coag_gain_a.txt" % coag_suffix, "w")
    outf_num_coag_gain_f = open("out/aging_%s_num_coag_gain_f.txt" % coag_suffix, "w")
    outf_num_coag_loss_a_a = open("out/aging_%s_num_coag_loss_a_a.txt" % coag_suffix, "w")
    outf_num_coag_loss_a_f = open("out/aging_%s_num_coag_loss_a_f.txt" % coag_suffix, "w")
    outf_num_coag_loss_f_a = open("out/aging_%s_num_coag_loss_f_a.txt" % coag_suffix, "w")
    outf_num_coag_loss_f_f = open("out/aging_%s_num_coag_loss_f_f.txt" % coag_suffix, "w")

    outf_mass_a = open("out/aging_%s_mass_a.txt" % coag_suffix, "w")
    outf_mass_f = open("out/aging_%s_mass_f.txt" % coag_suffix, "w")
    outf_mass_emit_a = open("out/aging_%s_mass_emit_a.txt" % coag_suffix, "w")
    outf_mass_emit_f = open("out/aging_%s_mass_emit_f.txt" % coag_suffix, "w")
    outf_mass_dilution_a = open("out/aging_%s_mass_dilution_a.txt" % coag_suffix, "w")
    outf_mass_dilution_f = open("out/aging_%s_mass_dilution_f.txt" % coag_suffix, "w")
    outf_mass_halving_a = open("out/aging_%s_mass_halving_a.txt" % coag_suffix, "w")
    outf_mass_halving_f = open("out/aging_%s_mass_halving_f.txt" % coag_suffix, "w")
    outf_mass_cond_a_a = open("out/aging_%s_mass_cond_a_a.txt" % coag_suffix, "w")
    outf_mass_cond_a_f = open("out/aging_%s_mass_cond_a_f.txt" % coag_suffix, "w")
    outf_mass_cond_f_a = open("out/aging_%s_mass_cond_f_a.txt" % coag_suffix, "w")
    outf_mass_cond_f_f = open("out/aging_%s_mass_cond_f_f.txt" % coag_suffix, "w")
    outf_mass_coag_gain_a = open("out/aging_%s_mass_coag_gain_a.txt" % coag_suffix, "w")
    outf_mass_coag_gain_f = open("out/aging_%s_mass_coag_gain_f.txt" % coag_suffix, "w")
    outf_mass_coag_loss_a_a = open("out/aging_%s_mass_coag_loss_a_a.txt" % coag_suffix, "w")
    outf_mass_coag_loss_a_f = open("out/aging_%s_mass_coag_loss_a_f.txt" % coag_suffix, "w")
    outf_mass_coag_loss_f_a = open("out/aging_%s_mass_coag_loss_f_a.txt" % coag_suffix, "w")
    outf_mass_coag_loss_f_f = open("out/aging_%s_mass_coag_loss_f_f.txt" % coag_suffix, "w")

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

        outf_height.write("%f %.20e" % (time, env_state.height))
        outf_comp_vol.write("%f %.20e" % (time, particles.comp_vol[0]))

        outf_num_a.write("%f " % time)
        outf_num_f.write("%f " % time)
        outf_num_emit_a.write("%f " % time)
        outf_num_emit_f.write("%f " % time)
        outf_num_dilution_a.write("%f " % time)
        outf_num_dilution_f.write("%f " % time)
        outf_num_halving_a.write("%f " % time)
        outf_num_halving_f.write("%f " % time)
        outf_num_cond_a_a.write("%f " % time)
        outf_num_cond_a_f.write("%f " % time)
        outf_num_cond_f_a.write("%f " % time)
        outf_num_cond_f_f.write("%f " % time)
        outf_num_coag_gain_a.write("%f " % time)
        outf_num_coag_gain_f.write("%f " % time)
        outf_num_coag_loss_a_a.write("%f " % time)
        outf_num_coag_loss_a_f.write("%f " % time)
        outf_num_coag_loss_f_a.write("%f " % time)
        outf_num_coag_loss_f_f.write("%f " % time)

        outf_mass_a.write("%f " % time)
        outf_mass_f.write("%f " % time)
        outf_mass_emit_a.write("%f " % time)
        outf_mass_emit_f.write("%f " % time)
        outf_mass_dilution_a.write("%f " % time)
        outf_mass_dilution_f.write("%f " % time)
        outf_mass_halving_a.write("%f " % time)
        outf_mass_halving_f.write("%f " % time)
        outf_mass_cond_a_a.write("%f " % time)
        outf_mass_cond_a_f.write("%f " % time)
        outf_mass_cond_f_a.write("%f " % time)
        outf_mass_cond_f_f.write("%f " % time)
        outf_mass_coag_gain_a.write("%f " % time)
        outf_mass_coag_gain_f.write("%f " % time)
        outf_mass_coag_loss_a_a.write("%f " % time)
        outf_mass_coag_loss_a_f.write("%f " % time)
        outf_mass_coag_loss_f_a.write("%f " % time)
        outf_mass_coag_loss_f_f.write("%f " % time)

        for ss_active in [0.001, 0.003, 0.006, 0.010]:
            num_a = 0
            num_f = 0
            num_emit_a = 0
            num_emit_f = 0
            num_dilution_a = 0
            num_dilution_f = 0
            num_halving_a = 0
            num_halving_f = 0
            num_cond_a_a = 0
            num_cond_a_f = 0
            num_cond_f_a = 0
            num_cond_f_f = 0
            num_coag_gain_a = 0
            num_coag_gain_f = 0
            num_coag_loss_a_a = 0
            num_coag_loss_a_f = 0
            num_coag_loss_f_a = 0
            num_coag_loss_f_f = 0

            mass_a = 0.0
            mass_f = 0.0
            mass_emit_a = 0.0
            mass_emit_f = 0.0
            mass_dilution_a = 0.0
            mass_dilution_f = 0.0
            mass_halving_a = 0.0
            mass_halving_f = 0.0
            mass_cond_a_a = 0.0
            mass_cond_a_f = 0.0
            mass_cond_f_a = 0.0
            mass_cond_f_f = 0.0
            mass_coag_gain_a = 0.0
            mass_coag_gain_f = 0.0
            mass_coag_loss_a_a = 0.0
            mass_coag_loss_a_f = 0.0
            mass_coag_loss_f_a = 0.0
            mass_coag_loss_f_f = 0.0

            if any(array(particles.aero_removed_action) == AERO_INFO_HALVED):
                halving_occured = True
            else:
                halving_occured = False

            # num, mass
            for i in range(particles.n_particles):
                if soot_mass[i] > 0.0:
                    if critical_ss[i] < ss_active:
                        num_a += 1
                        mass_a += soot_mass[i]
                    else:
                        num_f += 1
                        mass_f += soot_mass[i]

            if not first_time:
                removed_particles = {}
                for i in range(len(particles.aero_removed_id)):
                    removed_particles[particles.aero_removed_id[i]] = [particles.aero_removed_action[i],
                                                                       particles.aero_removed_other_id[i]]

                current_id = set()
                current_aged_id = set()
                current_fresh_id = set()
                for i in range(particles.n_particles):
                    current_id.add(particles.id[i])
                    if soot_mass[i] > 0.0:
                        if critical_ss[i] < ss_active:
                            current_aged_id.add(particles.id[i])
                        else:
                            current_fresh_id.add(particles.id[i])

                old_id = set()
                old_aged_id = set()
                old_fresh_id = set()
                for i in range(old_particles.n_particles):
                    old_id.add(old_particles.id[i])
                    if old_soot_mass[i] > 0.0:
                        if old_critical_ss[i] < ss_active:
                            old_aged_id.add(old_particles.id[i])
                        else:
                            old_fresh_id.add(old_particles.id[i])

                if not halving_occured:
                    if old_id - current_id != set(removed_particles.keys()):
                        print "old_id - current_id: ", old_id - current_id
                        print "removed_particles: ", set(removed_particles.keys())
                        print "xor: ", ((old_id - current_id) ^ set(removed_particles.keys()))
                        raise Exception("lost particle mismatch at t = %f" % time)

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
                for i in range(old_particles.n_particles):
                    if old_soot_mass[i] > 0.0:
                        id = old_particles.id[i]
                        if id in final_outcomes.keys():
                            [final_action, final_other_id] = final_outcomes[id]
                            if final_action == AERO_INFO_NONE:
                                raise Exception("found AERO_INFO_NONE at t = %f" % time)
                            elif final_action == AERO_INFO_DILUTION:
                                if old_critical_ss[i] < ss_active:
                                    num_dilution_a += 1
                                    mass_dilution_a += old_soot_mass[i]
                                else:
                                    num_dilution_f += 1
                                    mass_dilution_f += old_soot_mass[i]
                            elif final_action == AERO_INFO_COAG:
                                if final_other_id in current_aged_id:
                                    if old_critical_ss[i] < ss_active:
                                        num_coag_loss_a_a += 1
                                        mass_coag_loss_a_a += old_soot_mass[i]
                                    else:
                                        num_coag_loss_f_a += 1
                                        mass_coag_loss_f_a += old_soot_mass[i]
                                elif final_other_id in current_fresh_id:
                                    if old_critical_ss[i] < ss_active:
                                        num_coag_loss_a_f += 1
                                        mass_coag_loss_a_f += old_soot_mass[i]
                                    else:
                                        num_coag_loss_f_f += 1
                                        mass_coag_loss_f_f += old_soot_mass[i]
                                else:
                                    raise Exception("soot particle coagulated into a non-soot particle at t = %f" % time)
                            elif final_action == AERO_INFO_HALVED:
                                if old_critical_ss[i] < ss_active:
                                    num_halving_a += 1
                                    mass_halving_a += old_soot_mass[i]
                                else:
                                    num_halving_f += 1
                                    mass_halving_f += old_soot_mass[i]

                coag_id = set()
                for (id, [action, other_id]) in removed_particles.iteritems():
                    if action == AERO_INFO_COAG:
                        coag_id.add(id)
                        coag_id.add(other_id)

                if not (coag_id <= old_id):
                    raise Exception("coag_id not a subset of old_id at t = %f" % time)
                if not (coag_id <= set(final_outcomes.keys())):
                    raise Exception("coag_id not a subset of final_outcomes at t = %f" % time)

                # num_emit, num_cond, num_coag_gain
                for i in range(particles.n_particles):
                    id = particles.id[i]
                    if id in current_aged_id:
                        if id in coag_id:
                            num_coag_gain_a += 1
                            mass_coag_gain_a += soot_mass[i]
                        elif id in old_aged_id:
                            num_cond_a_a += 1
                            mass_cond_a_a += soot_mass[i]
                        elif id in old_fresh_id:
                            num_cond_f_a += 1
                            mass_cond_f_a += soot_mass[i]
                        elif id in old_id:
                            raise Exception("non-soot particle became soot particle without coagulation at t = %f" % time)
                        else:
                            num_emit_a += 1
                            mass_emit_a += soot_mass[i]
                    if id in current_fresh_id:
                        if id in coag_id:
                            num_coag_gain_f += 1
                            mass_coag_gain_f += soot_mass[i]
                        elif id in old_aged_id:
                            num_cond_a_f += 1
                            mass_cond_a_f += soot_mass[i]
                        elif id in old_fresh_id:
                            num_cond_f_f += 1
                            mass_cond_f_f += soot_mass[i]
                        elif id in old_id:
                            raise Exception("non-soot particle became soot particle without coagulation at t = %f" % time)
                        else:
                            num_emit_f += 1
                            mass_emit_f += soot_mass[i]
            
            outf_num_a.write("%d " % num_a)
            outf_num_f.write("%d " % num_f)
            outf_num_emit_a.write("%d " % num_emit_a)
            outf_num_emit_f.write("%d " % num_emit_f)
            outf_num_dilution_a.write("%d " % num_dilution_a)
            outf_num_dilution_f.write("%d " % num_dilution_f)
            outf_num_halving_a.write("%d " % num_halving_a)
            outf_num_halving_f.write("%d " % num_halving_f)
            outf_num_cond_a_a.write("%d " % num_cond_a_a)
            outf_num_cond_a_f.write("%d " % num_cond_a_f)
            outf_num_cond_f_a.write("%d " % num_cond_f_a)
            outf_num_cond_f_f.write("%d " % num_cond_f_f)
            outf_num_coag_gain_a.write("%d " % num_coag_gain_a)
            outf_num_coag_gain_f.write("%d " % num_coag_gain_f)
            outf_num_coag_loss_a_a.write("%d " % num_coag_loss_a_a)
            outf_num_coag_loss_a_f.write("%d " % num_coag_loss_a_f)
            outf_num_coag_loss_f_a.write("%d " % num_coag_loss_f_a)
            outf_num_coag_loss_f_f.write("%d " % num_coag_loss_f_f)
            
            outf_mass_a.write("%.20e " % mass_a)
            outf_mass_f.write("%.20e " % mass_f)
            outf_mass_emit_a.write("%.20e " % mass_emit_a)
            outf_mass_emit_f.write("%.20e " % mass_emit_f)
            outf_mass_dilution_a.write("%.20e " % mass_dilution_a)
            outf_mass_dilution_f.write("%.20e " % mass_dilution_f)
            outf_mass_halving_a.write("%.20e " % mass_halving_a)
            outf_mass_halving_f.write("%.20e " % mass_halving_f)
            outf_mass_cond_a_a.write("%.20e " % mass_cond_a_a)
            outf_mass_cond_a_f.write("%.20e " % mass_cond_a_f)
            outf_mass_cond_f_a.write("%.20e " % mass_cond_f_a)
            outf_mass_cond_f_f.write("%.20e " % mass_cond_f_f)
            outf_mass_coag_gain_a.write("%.20e " % mass_coag_gain_a)
            outf_mass_coag_gain_f.write("%.20e " % mass_coag_gain_f)
            outf_mass_coag_loss_a_a.write("%.20e " % mass_coag_loss_a_a)
            outf_mass_coag_loss_a_f.write("%.20e " % mass_coag_loss_a_f)
            outf_mass_coag_loss_f_a.write("%.20e " % mass_coag_loss_f_a)
            outf_mass_coag_loss_f_f.write("%.20e " % mass_coag_loss_f_f)

        outf_height.write("\n")
        outf_comp_vol.write("\n")

        outf_num_a.write("\n")
        outf_num_f.write("\n")
        outf_num_emit_a.write("\n")
        outf_num_emit_f.write("\n")
        outf_num_dilution_a.write("\n")
        outf_num_dilution_f.write("\n")
        outf_num_halving_a.write("\n")
        outf_num_halving_f.write("\n")
        outf_num_cond_a_a.write("\n")
        outf_num_cond_a_f.write("\n")
        outf_num_cond_f_a.write("\n")
        outf_num_cond_f_f.write("\n")
        outf_num_coag_gain_a.write("\n")
        outf_num_coag_gain_f.write("\n")
        outf_num_coag_loss_a_a.write("\n")
        outf_num_coag_loss_a_f.write("\n")
        outf_num_coag_loss_f_a.write("\n")
        outf_num_coag_loss_f_f.write("\n")

        outf_mass_a.write("\n")
        outf_mass_f.write("\n")
        outf_mass_emit_a.write("\n")
        outf_mass_emit_f.write("\n")
        outf_mass_dilution_a.write("\n")
        outf_mass_dilution_f.write("\n")
        outf_mass_halving_a.write("\n")
        outf_mass_halving_f.write("\n")
        outf_mass_cond_a_a.write("\n")
        outf_mass_cond_a_f.write("\n")
        outf_mass_cond_f_a.write("\n")
        outf_mass_cond_f_f.write("\n")
        outf_mass_coag_gain_a.write("\n")
        outf_mass_coag_gain_f.write("\n")
        outf_mass_coag_loss_a_a.write("\n")
        outf_mass_coag_loss_a_f.write("\n")
        outf_mass_coag_loss_f_a.write("\n")
        outf_mass_coag_loss_f_f.write("\n")

        old_particles = particles
        old_num_den = num_den
        old_critical_ss = critical_ss
        old_soot_mass = soot_mass
        first_time = False

    outf_height.close()
    outf_comp_vol.close()

    outf_num_a.close()
    outf_num_f.close()
    outf_num_emit_a.close()
    outf_num_emit_f.close()
    outf_num_dilution_a.close()
    outf_num_dilution_f.close()
    outf_num_halving_a.close()
    outf_num_halving_f.close()
    outf_num_cond_a_a.close()
    outf_num_cond_a_f.close()
    outf_num_cond_f_a.close()
    outf_num_cond_f_f.close()
    outf_num_coag_gain_a.close()
    outf_num_coag_gain_f.close()
    outf_num_coag_loss_a_a.close()
    outf_num_coag_loss_a_f.close()
    outf_num_coag_loss_f_a.close()
    outf_num_coag_loss_f_f.close()
    
    outf_mass_a.close()
    outf_mass_f.close()
    outf_mass_emit_a.close()
    outf_mass_emit_f.close()
    outf_mass_dilution_a.close()
    outf_mass_dilution_f.close()
    outf_mass_halving_a.close()
    outf_mass_halving_f.close()
    outf_mass_cond_a_a.close()
    outf_mass_cond_a_f.close()
    outf_mass_cond_f_a.close()
    outf_mass_cond_f_f.close()
    outf_mass_coag_gain_a.close()
    outf_mass_coag_gain_f.close()
    outf_mass_coag_loss_a_a.close()
    outf_mass_coag_loss_a_f.close()
    outf_mass_coag_loss_f_a.close()
    outf_mass_coag_loss_f_f.close()
