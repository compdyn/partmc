#!/usr/bin/env python

import sys
sys.path.append('../../../../tool')
import partmc
import scipy
import scipy.io
ncf = scipy.io.netcdf_file('condense_0001_00000061.nc')
p = partmc.aero_particle_array_t(ncf)
n_part = len(p.num_concs)
n_ccn = len(p.num_concs[p.diameters() > 2e-6])
conc_part = p.num_concs.sum()
conc_ccn = (p.num_concs[p.diameters() > 2e-6]).sum()
ccn_cn = conc_ccn / conc_part

print "n_part = ", n_part
print "n_ccn = ", n_ccn
print "conc_part = ", conc_part
print "conc_ccn = ", conc_ccn
print "ccn_cn = ", ccn_cn
