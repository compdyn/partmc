#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math, re
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *
import numpy
sys.path.append(".")
from fig_helper import *

data_prefix = "data/particles"

particle_ids = [p["id"] for p in show_particles]
particle_histories = read_history(lambda ncf:
                                  aero_particle_array_t(ncf, include_ids = particle_ids),
                                  netcdf_dir_wc, netcdf_pattern_wc, print_progress = True)
particle_histories.sort()

for id in particle_ids:
    print id
    aero_data = read_any(aero_data_t, netcdf_dir_wc, netcdf_pattern_wc)
    particle_history = [(t, p) for (t, p) in particle_histories
                        if id in p.id]
    if len(particle_history) == 0:
        print "WARNING: particle ID %d not found" % id
        continue
    data = zeros([len(particle_history),
                  len(aero_data.name) + 2])
    for (i, (time, particles)) in enumerate(particle_history):
        index = nonzero(particles.id == id)[0][0]
        data[i, 0] = time
        data[i, 1] = particles.comp_vol[index]
        for s in range(len(aero_data.name)):
            data[i, s+2] = particles.masses[s, index]
    filename = "%s_%d.txt" % (data_prefix, id)
    savetxt(filename, data)
