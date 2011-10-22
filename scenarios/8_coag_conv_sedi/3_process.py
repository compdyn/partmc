#!/usr/bin/env python

import sys, os, numpy, scipy
sys.path.append("../../tool")
import partmc
import config

for run in config.all_runs():
    dirname = os.path.join(config.run_dirname, run["name"])
    print dirname

    sect_filename = os.path.join(dirname, "sect_00000002.nc")
    ncf = scipy.io.netcdf.netcdf_file(sect_filename, 'r')
    diam_grid = partmc.log_grid()
    diam_grid.load_ncf_diam(ncf)
    aero_binned = partmc.aero_binned_t(ncf)
    ncf.close()

    part_filenames = partmc.get_filename_list(dirname, r"part_[0-9]*_00000002\.nc")
    num_err = numpy.zeros(len(part_filenames))
    mass_err = numpy.zeros(len(part_filenames))
    for (i, part_filename) in enumerate(part_filenames):
        ncf = scipy.io.netcdf.netcdf_file(part_filename, 'r')
        aero_particle_array = partmc.aero_particle_array_t(ncf)
        ncf.close()

        diameters = aero_particle_array.diameters()
        masses = aero_particle_array.masses()
        comp_vols = aero_particle_array.comp_vols

        num_dist = partmc.histogram_1d(diameters, diam_grid, weights=1 / comp_vols)
        mass_dist = partmc.histogram_1d(diameters, diam_grid, weights=masses / comp_vols)

        num_err[i] = numpy.linalg.norm(num_dist - aero_binned.num_conc)
        mass_err[i] = numpy.linalg.norm(mass_dist - aero_binned.mass_conc())

    num_err_mean = num_err.mean()
    num_err_ci = 2 * num_err.std() / numpy.sqrt(len(num_err))
    print "num_err", num_err
    print "num_mean", num_err.mean()
    print "num_std", num_err.std()
    print "num_ci", num_err_ci
    mass_err_mean = mass_err.mean()
    mass_err_ci = 2 * mass_err.std() / numpy.sqrt(len(mass_err))
    print "mass_err", mass_err
    print "mass_mean", mass_err.mean()
    print "mass_std", mass_err.std()
    print "mass_ci", mass_err_ci

    stats = numpy.zeros(4)
    stats[0] = num_err_mean
    stats[1] = num_err_ci
    stats[2] = mass_err_mean
    stats[3] = mass_err_ci

    stats_filename = os.path.join(dirname, "stats.txt")
    numpy.savetxt(stats_filename, stats)
