#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math, re
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
import numpy
sys.path.append(".")
from fig_helper import *

out_prefix = "out/aging_aero_particles"

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)

plot_data = [[] for p in show_particles]
for [time, filename, key] in time_filename_list:
    print filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    for (pi, p) in enumerate(show_particles):
        i = (abs(particles.id - p["id"])).argmin()
        if particles.id[i] == p["id"]:
            print p["id"]
            #if filename == "../urban_plume/out/urban_plume_wc_state_0001_00000372.nc":
            #    print particles.mass
            n_spec = size(particles.masses, 0)
            masses = [particles.masses[s,i] for s in range(n_spec)]
            plot_data[pi].append([time, masses])

for (pi, p) in enumerate(show_particles):
    print p["id"]
    data = zeros([len(plot_data[pi]), n_spec + 1])
    for (i_time, [time, masses]) in enumerate(plot_data[pi]):
        data[i_time, 0] = time
        for i_spec in range(n_spec):
            data[i_time, 1 + i_spec] = masses[i_spec]
    savetxt("%s_%d.txt" % (out_prefix, pi), data)
