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

data_prefix = "data/particles"

particle_ids = [p["id"] for p in show_particles]
for id in particle_ids:
    aero_data = read_any(aero_data_t, netcdf_dir_wc, netcdf_pattern_wc)
    particle_history = read_history(lambda ncf:
                                    aero_particle_array_t(ncf, include_ids = [id]),
                                    netcdf_dir_wc, netcdf_pattern_wc)
    particle_history = [(t, p) for (t, p) in particle_history
                        if size(p.id) > 0]
    if len(particle_history) == 0:
        print "WARNING: particle ID %d not found" % id
    data = array([len(particle_history),
                  len(aero_data.name) + 2])
    for (i, (time, particles)) in enumerate(particle_history):
        data[i, 0] = time
        data[i, 1] = particles.comp_vol[0]
        for s in range(len(aero_data.name)):
            data[i, s+2] = particles.masses[s, 0]
    filename = "%s_%d.txt" % (data_prefix, id)
    savetxt(filename, data)
