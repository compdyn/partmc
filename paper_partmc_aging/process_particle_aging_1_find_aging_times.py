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

    time_emitted = {}
    mass_emitted = {}
    time_entered_bins = [{} for i in range(n_level_bin + 2)]

    first_time = True
    for [time, filename, key] in time_filename_list:
        #DEBUG
        if time > 121:
            sys.exit(0)
        #DEBUG
        print time, filename
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        particles.id = [int(i) for i in particles.id]
        env_state = env_state_t(ncf)
        ncf.close()
        if first_time:
            aero_data = particles.aero_data
            n_spec = len(aero_data.name)
        num_den = 1.0 / array(particles.comp_vol)
        total_num_den = num_den.sum()
        soot_mass = particles.mass(include = ["BC"])
        critical_ss = particles.kappa_rh(env_state, const) - 1.0
        ss_bin = ss_active_axis.find_clipped_outer(critical_ss)
        id_set = set([particles.id[i] for i in range(particles.n_particles)])
        particle_index_by_id = dict([[particles.id[i], i] for i in range(particles.n_particles)])

        for i in range(particles.n_particles):
            id = particles.id[i]
            if id not in time_emitted:
                time_emitted[id] = time
                mass_emitted[id] = [particles.masses[s, i] for s in range(size(particles.masses, 0))]
            bin = ss_bin[i]
            if id not in time_entered_bins[bin]:
                time_entered_bins[bin][id] = time

    max_id = max(time_emitted.keys())

    filename = os.path.join(aging_data_dir,
                            "particle_aging_%s_%%s.txt" % coag_suffix)

    time_emitted_array = zeros(max_id)
    mass_emitted_array = zeros([max_id, n_spec])
    time_emitted_array = -1.0
    mass_emitted_array = -1.0
    for (id, time) in time_emitted.iteritems():
        time_emitted_array[id] = time
        for s in range(n_spec):
            mass_emitted_array[id, s] = mass_emitted[id][s]
    savetxt(filename % "time_emitted", time_emitted_array)
    savetxt(filename % "mass_emitted", mass_emitted_array)

    for bin in range(n_level_bin + 2):
        print "bin", bin
        filename_bin = filename % ("%08d_%%s" % bin)

        time_entered_bins_array = zeros(max_id)
        time_aging_array = zeros(max_id)

        time_entered_bins_array = -1.0
        time_aging_array = -1.0

        for (id, time) in time_entered_bins[bin]:
            time_entered_bins_array[id] = time
            if time_emitted_array[id] < 0.0:
                raise Exception("ID %d entered bin %d at time %f but was never emitted" % id, bin, time)
            time_aging_array[id] = time - time_emitted_array[id]

        savetxt(filename_bin % "time_entered_bins", time_entered_bins_array)
        savetxt(filename_bin % "time_aging", time_aging_array)
