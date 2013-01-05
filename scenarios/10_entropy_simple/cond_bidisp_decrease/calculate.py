#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
sys.path.append("../../../tool")
import partmc

ncf = scipy.io.netcdf_file("out/urban_plume_0001_00000003.nc")

particles = partmc.aero_particle_array_t(ncf)

dry_diameters = particles.dry_diameters()
diameters = particles.diameters()

bc = particles.masses(include = ["BC"])
so4 = particles.masses(include = ["SO4"])
nh4 = particles.masses(include = ["NH4"])
no3 = particles.masses(include = ["NO3"])
h2o = particles.masses(include = ["H2O"])

contains_no_so4 = (so4 == 0)
number_mode_1 = sum(contains_no_so4)

contains_so4 = (so4 > 0)
number_mode_2 = sum(contains_so4)

num_concs = particles.num_concs
bc_total = sum(bc*num_concs)
so4_total = sum(so4*num_concs)
nh4_total = sum(nh4*num_concs)
no3_total = sum(no3*num_concs)
h2o_total = sum(h2o*num_concs)

num_concs_mode_1 = sum(particles.num_concs[contains_no_so4])
bc_mode_1 = sum(bc[contains_no_so4]*num_concs[contains_no_so4])
so4_mode_1 = sum(so4[contains_no_so4]*num_concs[contains_no_so4])
nh4_mode_1 = sum(nh4[contains_no_so4]*num_concs[contains_no_so4])
no3_mode_1 = sum(no3[contains_no_so4]*num_concs[contains_no_so4])
h2o_mode_1 = sum(h2o[contains_no_so4]*num_concs[contains_no_so4])

num_concs_mode_2 = sum(particles.num_concs[contains_so4])
bc_mode_2 = sum(bc[contains_so4]*num_concs[contains_so4])
so4_mode_2 = sum(so4[contains_so4]*num_concs[contains_so4])
nh4_mode_2 = sum(nh4[contains_so4]*num_concs[contains_so4])
no3_mode_2 = sum(no3[contains_so4]*num_concs[contains_so4])
h2o_mode_2 = sum(h2o[contains_so4]*num_concs[contains_so4])

mass_mode_1 = bc_mode_1+so4_mode_1+nh4_mode_1+no3_mode_1+h2o_mode_1
mass_mode_2 = bc_mode_2+so4_mode_2+nh4_mode_2+no3_mode_2+h2o_mode_2

diameters_1 = diameters[contains_no_so4]
diameters_2 = diameters[contains_so4]

print "mode 1 (pure BC}"
print "num_conc", num_concs_mode_1
print "total mass ", mass_mode_1
print "bc ", bc_mode_1, bc_mode_1/mass_mode_1
print "so4 ", so4_mode_1, so4_mode_1/mass_mode_1
print "nh4 ", nh4_mode_1, nh4_mode_1/mass_mode_1
print "no3 ", no3_mode_1, no3_mode_1/mass_mode_1
print "h2o ", h2o_mode_1, h2o_mode_1/mass_mode_1
print "total ", bc_mode_1/mass_mode_1+so4_mode_1/mass_mode_1+nh4_mode_1/mass_mode_1+ no3_mode_1/mass_mode_1+h2o_mode_1/mass_mode_1
print "diameter ", diameters_1[1]
print
print "mode 2 (mixed BC}"
print "num_conc", num_concs_mode_2
print "total mass ", mass_mode_2
print "bc ", bc_mode_2, bc_mode_2/mass_mode_2
print "so4 ", so4_mode_2, so4_mode_2/mass_mode_2
print "nh4 ", nh4_mode_2, nh4_mode_2/mass_mode_2
print "no3 ", no3_mode_2, no3_mode_2/mass_mode_2
print "h2o ", h2o_mode_2, h2o_mode_2/mass_mode_2
print "total ", bc_mode_2/mass_mode_2+so4_mode_2/mass_mode_2+nh4_mode_2/mass_mode_2+ no3_mode_2/mass_mode_2+h2o_mode_2/mass_mode_2
print "diameter ", diameters_2[1] 





