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

bin = level_mid + 1

for time_type in ["emission", "aging"]:
    for coag in [True]:
        if coag:
            env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
            aero_data = read_any(aero_data_t, netcdf_dir_wc, netcdf_pattern_wc)
            coag_suffix = "wc"
        else:
            env_state = read_any(env_state_t, netcdf_dir_nc, netcdf_pattern_nc)
            aero_data = read_any(aero_data_t, netcdf_dir_nc, netcdf_pattern_nc)
            coag_suffix = "nc"
        print time_type, coag_suffix

        filename = os.path.join(aging_data_dir,
                                "particle_aging_%s_%%s.txt" % coag_suffix)
        filename_bin = filename % ("%08d_%%s" % bin)
        time_emitted_array = loadtxt(filename % "time_emitted")
        mass_emitted_array = loadtxt(filename % "mass_emitted")
        time_aging_array = loadtxt(filename_bin % "time_aging")

        n_particles = 0
        for id in range(size(time_emitted_array)):
            if time_emitted_array[id] >= 0.0:
                n_particles += 1
        particles = aero_particle_array_t(n_particles = n_particles,
                                          aero_data = aero_data)
        i = 0
        particle_index_by_id = {}
        for id in range(size(time_emitted_array)):
            if time_emitted_array[id] >= 0.0:
                particle_index_by_id[id] = i
                particles.id[i] = id
                for s in range(len(aero_data.name)):
                    particles.masses[s,i] = mass_emitted_array[id,s]
                i += 1

        diameter = particles.dry_diameter() * 1e6
        bc_mass = particles.mass(include = ["BC"])

        num_den_array = numpy.zeros([diameter_axis.n_bin, aging_time_axis.n_bin])

        aging_time = zeros_like(diameter)
        aging_time[:] = -1.0
        for id in range(size(time_emitted_array)):
            if time_aging_array[id] >= 0.0:
                i = particle_index_by_id[id]
                aging_time[i] = time_aging_array[id]
        aging_time /= 3600.0

        if time_type == "emission":
            z_data = zeros_like(diameter)
            for i in range(size(diameter)):
                id = particles.id[i]
                z_data[i] = time_emitted_array[id] / 3600.0
        elif time_type == "aging":
            z_data = zeros_like(diameter)
            for i in range(size(diameter)):
                id = particles.id[i]
                aging_time_point = time_emitted_array[id] + time_aging_array[id]
                z_data[i] = aging_time_point / 3600.0
        else:
            raise Exception("unknown time_type: %s" % time_type)

        n_use = 0
        for i in range(size(diameter)):
            if bc_mass[i] > 0.0:
                n_use += 1

        use_data = zeros([n_use, 3])

        use_i = 0
        for i in range(size(diameter)):
            if bc_mass[i] > 0.0:
                use_data[use_i, 0] = diameter[i]
                use_data[use_i, 1] = aging_time[i]
                use_data[use_i, 2] = z_data[i]
                use_i += 1

        filename = os.path.join(aging_data_dir,
                                "particle_aging_by_time_%s_%s_plot_data_%08d.txt" % (time_type, coag_suffix, bin))

        savetxt(filename, use_data)
