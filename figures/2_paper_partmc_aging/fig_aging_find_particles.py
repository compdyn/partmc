#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math, random
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append(".")
from fig_helper import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

const = load_constants("../src/constants.f90")

time_hours = [12, 15, 21, 24]

def read_particles(filename):
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()
    return (particles, env_state)

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)

# particles that never coagulated
(_, filename, _) = time_filename_list[-1]
(particles, env_state) = read_particles(filename)
no_coag_ids = set()
for i in range(particles.n_particles):
    if particles.n_orig_part[i] == 1:
        no_coag_ids.add(particles.id[i])

def find_aged_or_fresh_ids(particles, env_state, aged):
    critical_ss = particles.kappa_rh(env_state, const) - 1.0
    ids = set()
    for i in range(particles.n_particles):
        if (aged and critical_ss[i] < level_mid_value) \
                or (not aged and critical_ss[i] > level_mid_value):
            if particles.id[i] in no_coag_ids:
                ids.add(particles.id[i])
    return ids

time_hours = [15, 18, 21, 24]
particles = {}
env_states = {}
for time_hour in time_hours:
    filename = file_filename_at_time(time_filename_list, time_hour * 2400.0)
    (particles[time_hour], env_states[time_hour]) = read_particles(filename)

aged_15_ids = find_aged_or_fresh_ids(particles[15], env_states[15], True)
fresh_18_ids = find_aged_or_fresh_ids(particles[18], env_states[18], False)
aged_21_ids = find_aged_or_fresh_ids(particles[21], env_states[21], True)
fresh_24_ids = find_aged_or_fresh_ids(particles[24], env_states[24], False)

deaged_15_18_ids = aged_15_ids & fresh_18_ids
deaged_21_24_ids = aged_21_ids & fresh_24_ids
deaged_both_ids = deaged_15_18_ids & deaged_21_24_ids

print "de-aged 15 -> 18 hours:"
diameter_15 = particles[15].dry_diameter() * 1e6
create_time_15 = particles[15].least_create_time / 3600.0
data = []
for id in deaged_15_18_ids:
    i = (abs(particles[15].id - id)).argmin()
    data.append((id, diameter_15[i], create_time_15[i]))
data.sort(key = lambda x: x[2])
for (id, diameter, create_time) in data:
    print "id: %d --- diameter: %f um --- create time: %f h" \
        % (id, diameter, create_time)
    
print "de-aged 21 -> 24 hours:"
diameter_21 = particles[21].dry_diameter() * 1e6
create_time_21 = particles[21].least_create_time / 3600.0
data = []
for id in deaged_21_24_ids:
    i = (abs(particles[21].id - id)).argmin()
    data.append((id, diameter_21[i], create_time_21[i]))
data.sort(key = lambda x: x[2])
for (id, diameter, create_time) in data:
    print "id: %d --- diameter: %f um --- create time: %f h" \
        % (id, diameter, create_time)

print "de-aged both:"
print deaged_both_ids
