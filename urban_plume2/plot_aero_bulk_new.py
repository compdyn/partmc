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

aero_species = [["", ["NO3"]],
                ["", ["NH4"]],
                ["", ["OC"]],
                ["", ["H2O"]],
                ["", ["SO4"]],
                ["", ["BC"]],
                ["SOA", ["ARO1", "ARO2", "ALK1", "OLE1"]],
                ]

netcdf_dir = "out"
netcdf_pattern = r"urban_plume_state_0001_([0-9]{8})\.nc"

filenames = os.listdir(netcdf_dir)
num_data = []
mass_data = []
mass_species_data = [[] for x in aero_species]
max_time_min = 0.0
filename_re = re.compile(netcdf_pattern)
for filename in filenames:
    if filename_re.search(filename):
        netcdf_filename = os.path.join(netcdf_dir, filename)
        print netcdf_filename
        ncf = NetCDFFile(netcdf_filename)
        env_state = env_state_t(ncf)
        time_min = env_state.elapsed_time / 60
        max_time_min = max(max_time_min, time_min)
        gas_state = gas_state_t(ncf)
        particle_array = aero_particle_array_t(ncf)
        
        num_den = (1.0 / particle_array.comp_vol).sum()
        if num_den > 0.0:
            num_data.append([time_min, num_den])
        
        total_masses = particle_array.mass()
        total_mass_den = (total_masses / particle_array.comp_vol).sum()
        if total_mass_den > 0.0:
            mass_data.append([time_min, total_mass_den])
        
        for i in range(len(aero_species)):
            species_masses = particle_array.mass(include = aero_species[i][1])
            species_mass_den = (species_masses / particle_array.comp_vol).sum()
            if species_mass_den > 0.0:
                mass_species_data[i].append([time_min, species_mass_den])
            
        ncf.close()
num_data.sort()
mass_data.sort()
for i in range(len(aero_species)):
    mass_species_data[i].sort

env_state = read_any(env_state_t, netcdf_dir, netcdf_pattern)
start_time_of_day_min = env_state.start_time_of_day / 60

######################################################################

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.0,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "time (LST)",
			  painter = grid_painter),
    y = graph.axis.log(title = r"number density ($1/m^3$)",
                       painter = grid_painter))

g.plot(graph.data.points(num_data, x = 1, y = 2),
       styles = [graph.style.line(lineattrs = [color.rgb.black,
                                               style.linewidth.THick])])

g.writePDFfile("out/aero_bulk_num.pdf")

######################################################################

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.0,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "time (LST)",
			  painter = grid_painter),
    y = graph.axis.log(title = r"mass density ($m^3/m^3$)",
                       painter = grid_painter))

g.plot(graph.data.points(mass_data, x = 1, y = 2),
       styles = [graph.style.line(lineattrs = [color.rgb.black,
                                               style.linewidth.THick])])

g.writePDFfile("out/aero_bulk_mass.pdf")

######################################################################

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.0,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "time (LST)",
			  painter = grid_painter),
    y = graph.axis.log(title = r"mass density ($m^3/m^3$)",
                       painter = grid_painter),
    key = graph.key.key(pos = "br"))

for i in range(len(aero_species)):
    if aero_species[i][0] == "":
        title = tex_species(aero_species[i][1][0])
        print_title = aero_species[i][1][0]
    else:
        title = aero_species[i][0]
        print_title = aero_species[i][0]
    if mass_species_data[i] != []:
        g.plot(graph.data.points(mass_species_data[i],
                                 x = 1, y = 2, title = title),
               styles = [graph.style.line(lineattrs = [color_list[i],
                                                       style.linewidth.THick])])
    else:
        print "warning: only zeros for species %s" % print_title

g.writePDFfile("out/aero_bulk_by_species.pdf")

######################################################################
