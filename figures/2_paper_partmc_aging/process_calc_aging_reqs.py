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

aero_data = read_any(aero_data_t, netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
print "env_state.time = ", env_state.elapsed_time
env_state.temperature = 273.0
print "env_state.temperature = ", env_state.temperature

particles = aero_particle_array_t(n_particles = 1, aero_data = aero_data)

def set_bc_mass_by_diameter(desired_d):
    particles.masses[:,:] = 0.0
    def f(bc_mass):
        particles.masses[aero_data.name.index("BC")] = bc_mass
        diameter = particles.diameter()
        error = diameter[0] - desired_d
        return error
    bc_mass = scipy.optimize.brentq(f, 0.0, 1.0, xtol = 1e-30)
    particles.masses[aero_data.name.index("BC")] = bc_mass

def set_nh4no3_by_css(desired_css):
    def f(nh4_mass):
        nh4_moles = nh4_mass / aero_data.molec_weight[aero_data.name.index("NH4")]
        no3_moles = nh4_moles
        no3_mass = no3_moles * aero_data.molec_weight[aero_data.name.index("NO3")]
        particles.masses[aero_data.name.index("NH4")] = nh4_mass
        particles.masses[aero_data.name.index("NO3")] = no3_mass
        critical_ss = particles.kappa_rh(env_state, const) - 1.0
        error = critical_ss[0] - desired_css
        return error
    bc_mass = particles.masses[aero_data.name.index("BC")]
    nh4_mass = scipy.optimize.brentq(f, bc_mass * 1e-5, 1.0, xtol = 1e-30)
    nh4_moles = nh4_mass / aero_data.molec_weight[aero_data.name.index("NH4")]
    no3_moles = nh4_moles
    no3_mass = no3_moles * aero_data.molec_weight[aero_data.name.index("NO3")]
    particles.masses[aero_data.name.index("NH4")] = nh4_mass
    particles.masses[aero_data.name.index("NO3")] = no3_mass

def print_info(with_css = False):
    print "diameter = ", particles.diameter()[0]
    print "BC mass = ", particles.masses[aero_data.name.index("BC"), 0]
    print "NH4 mass = ", particles.masses[aero_data.name.index("NH4"), 0]
    print "NO3 mass = ", particles.masses[aero_data.name.index("NO3"), 0]
    print "NH4NO3 mass = %.3f fg" % ((particles.masses[aero_data.name.index("NH4"), 0]
                                     + particles.masses[aero_data.name.index("NO3"), 0]) * 1e18)
    print "kappa = ", particles.solute_kappa()[0]
    if with_css:
        print "critical supersaturation = ", particles.kappa_rh(env_state, const)[0] - 1.0

print "**************************************"
set_bc_mass_by_diameter(0.02e-6)
print_info()
print "---------------------------"
set_nh4no3_by_css(0.3e-2)
print_info(with_css = True)

print "**************************************"
set_bc_mass_by_diameter(0.2e-6)
print_info()
print "---------------------------"
set_nh4no3_by_css(0.3e-2)
print_info(with_css = True)
