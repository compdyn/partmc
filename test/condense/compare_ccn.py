#!/usr/bin/env python

import sys
sys.path.append('../../../../tool')
import partmc
import scipy.io
import numpy

ncf = scipy.io.netcdf_file('condense_0001_00000601.nc')
p = partmc.aero_particle_array_t(ncf)
ncf.close()
p.num_concs = 1 / p.comp_vols
n_part = len(p.num_concs)
n_ccn = len(p.num_concs[p.diameters() > 2e-6])
conc_part = p.num_concs.sum()
conc_ccn = (p.num_concs[p.diameters() > 2e-6]).sum()
ccn_cn = conc_ccn / conc_part

n_sel = 20
ncf = scipy.io.netcdf_file('condense_0001_00000001.nc')
p = partmc.aero_particle_array_t(ncf)
ncf.close()
d = p.diameters()
desired_diams = numpy.logspace(numpy.log10(d.min()), numpy.log10(d.max()), n_sel)
ids = numpy.zeros(n_sel)
for i in range(n_sel):
    idx = numpy.abs(d - desired_diams[i]).argmin()
    ids[i] = p.ids[idx]

rh_max = 0
n_time = 601
data = numpy.zeros((n_time, n_sel + 1))
for i in range(n_time):
    print i
    filename = 'condense_0001_%08d.nc' % (i + 1)
    ncf = scipy.io.netcdf_file(filename)
    env_state = partmc.env_state_t(ncf)
    p = partmc.aero_particle_array_t(ncf)
    ncf.close()
    rh_max = max(rh_max, env_state.relative_humidity)
    d = p.diameters()
    data[i, 0] = env_state.elapsed_time
    for j in range(n_sel):
        idx = numpy.abs(p.ids - ids[j]).argmin()
        data[i, j + 1] = d[idx]
print ids

numpy.savetxt("diams.txt", data)

print "n_part = ", n_part
print "n_ccn = ", n_ccn
print "conc_part = ", conc_part
print "conc_ccn = ", conc_ccn
print "ccn_cn = ", ccn_cn
print "rh_max = ", rh_max
