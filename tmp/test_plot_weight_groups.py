#!/usr/bin/env python

import os, sys
import numpy
import scipy.io
sys.path.append("../../../tool")
import partmc

def process_file(filename, prefix):
    ncf = scipy.io.netcdf.netcdf_file(filename)
    particles = partmc.aero_particle_array_t(ncf)
    n_groups = ncf.dimensions["aero_weight"]
    comp_vols = ncf.variables["weight_comp_vol"][()]
    exponents = ncf.variables["weight_exponent"][()]
    weight_types = ncf.variables["weight_type"][()]
    ncf.close()

    x_axis = partmc.log_grid(min=1e-7,max=1e-2,n_bin=100)
    x_centers = x_axis.centers()
    num = numpy.zeros([len(x_centers), n_groups + 4])
    mass = numpy.zeros([len(x_centers), n_groups + 4])
    V_loc_sum = numpy.zeros(len(x_centers))

    diameters = particles.diameters()
    masses = particles.masses()
    num[:,0] = x_centers
    mass[:,0] = x_centers
    num[:,1] = partmc.histogram_1d(diameters, x_axis,
                                   weights = 1 / particles.comp_vols) / numpy.log(10)
    mass[:,1] = partmc.histogram_1d(diameters, x_axis,
                                    weights = masses / particles.comp_vols) / numpy.log(10)
    for i_group in range(1, n_groups + 1):
        this_group = (particles.weight_groups == i_group)
        num[:, i_group + 3] = partmc.histogram_1d(diameters[this_group],
                                              x_axis, weights = 1 / particles.comp_vols[this_group]) / numpy.log(10)
        mass[:, i_group + 3] = partmc.histogram_1d(diameters[this_group], x_axis,
                                               weights = masses[this_group] / particles.comp_vols[this_group]) / numpy.log(10)
        bin_weights = x_centers**exponents[i_group - 1]
        V_loc = comp_vols[i_group - 1] / bin_weights
        V_loc_sum += V_loc
        num[:,2] += V_loc * num[:, i_group + 3]
        mass[:,2] += V_loc * mass[:, i_group + 3]

        num[:,3] += num[:, i_group + 3]
        mass[:,3] += mass[:, i_group + 3]

    num[:,2] *= n_groups / V_loc_sum
    mass[:,2] *= n_groups / V_loc_sum

    numpy.savetxt(prefix + "num.txt", num)
    numpy.savetxt(prefix + "mass.txt", mass)

process_file("out/sedi_part_0001_00000001.nc", "out/1_")
process_file("out/sedi_part_0001_00000002.nc", "out/2_")
process_file("out/sedi_part_0001_00000003.nc", "out/3_")
