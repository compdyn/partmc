#!/usr/bin/env python

import sys, os, numpy, scipy
sys.path.append("../../tool")
import partmc
import config

sect_filename = os.path.join(config.run_dirname, "1k_flat_source", "sect_00000002.nc")
ncf = scipy.io.netcdf.netcdf_file(sect_filename, 'r')
diam_grid = partmc.log_grid()
diam_grid.load_ncf_diam(ncf)
aero_binned_array = partmc.aero_binned_array_t(ncf)
ncf.close()

sect_filename = os.path.join(config.run_dirname, "1k_flat", "sect_00000002.nc")
ncf = scipy.io.netcdf.netcdf_file(sect_filename, 'r')
aero_binned_single = partmc.aero_binned_t(ncf)
ncf.close()

for run in config.all_runs():
    dirname = os.path.join(config.run_dirname, run["name"])
    print dirname

    part_filenames = partmc.get_filename_list(dirname, r"part_[0-9]{4}_00000002\.nc")
    num_1_err = numpy.zeros(len(part_filenames))
    num_2_err = numpy.zeros(len(part_filenames))
    mass_1_err = numpy.zeros(len(part_filenames))
    mass_2_err = numpy.zeros(len(part_filenames))
    num_1_dists = numpy.zeros((len(part_filenames), diam_grid.n_bin))
    num_2_dists = numpy.zeros((len(part_filenames), diam_grid.n_bin))
    mass_1_dists = numpy.zeros((len(part_filenames), diam_grid.n_bin))
    mass_2_dists = numpy.zeros((len(part_filenames), diam_grid.n_bin))
    num_dists = numpy.zeros((len(part_filenames), diam_grid.n_bin))
    mass_dists = numpy.zeros((len(part_filenames), diam_grid.n_bin))
    for (i, part_filename) in enumerate(part_filenames):
        ncf = scipy.io.netcdf.netcdf_file(part_filename, 'r')
        aero_particle_array = partmc.aero_particle_array_t(ncf)
        ncf.close()

        diameters = aero_particle_array.diameters()
        masses = aero_particle_array.masses()
        masses2 = aero_particle_array.masses(include=["S2"])
        comp_vols = aero_particle_array.comp_vols

        species1 = (masses2 == 0)
        species2 = (masses2 > 0)

        num_dist = partmc.histogram_1d(diameters, diam_grid, weights=1 / comp_vols)
        mass_dist = partmc.histogram_1d(diameters, diam_grid, weights=masses / comp_vols)

        num_1_dist = partmc.histogram_1d(diameters[species1], diam_grid, weights=1 / comp_vols[species1])
        num_2_dist = partmc.histogram_1d(diameters[species2], diam_grid, weights=1 / comp_vols[species2])
        mass_1_dist = partmc.histogram_1d(diameters[species1], diam_grid, weights=masses[species1] / comp_vols[species1])
        mass_2_dist = partmc.histogram_1d(diameters[species2], diam_grid, weights=masses[species2] / comp_vols[species2])

        num_1_err[i] = numpy.linalg.norm(num_1_dist - aero_binned_array.num_conc[0,:])
        num_2_err[i] = numpy.linalg.norm(num_2_dist - aero_binned_array.num_conc[1,:])
        mass_1_err[i] = numpy.linalg.norm(mass_1_dist - aero_binned_array.mass_conc()[0,:])
        mass_2_err[i] = numpy.linalg.norm(mass_2_dist - aero_binned_array.mass_conc()[1,:])

        num_dists[i,:] = num_dist
        mass_dists[i,:] = mass_dist

        num_1_dists[i,:] = num_1_dist
        num_2_dists[i,:] = num_2_dist
        mass_1_dists[i,:] = mass_1_dist
        mass_2_dists[i,:] = mass_2_dist

    num_1_err_mean = num_1_err.mean()
    num_1_err_ci = 2 * num_1_err.std() / numpy.sqrt(len(num_1_err))
    mass_1_err_mean = mass_1_err.mean()
    mass_1_err_ci = 2 * mass_1_err.std() / numpy.sqrt(len(mass_1_err))

    num_2_err_mean = num_2_err.mean()
    num_2_err_ci = 2 * num_2_err.std() / numpy.sqrt(len(num_2_err))
    mass_2_err_mean = mass_2_err.mean()
    mass_2_err_ci = 2 * mass_2_err.std() / numpy.sqrt(len(mass_2_err))

    num_1_dist_mean = num_1_dists.mean(axis=0)
    num_1_dist_ci = 2 * num_1_dists.std(axis=0)
    num_2_dist_mean = num_2_dists.mean(axis=0)
    num_2_dist_ci = 2 * num_2_dists.std(axis=0)
    mass_1_dist_mean = mass_1_dists.mean(axis=0)
    mass_1_dist_ci = 2 * mass_1_dists.std(axis=0)
    mass_2_dist_mean = mass_2_dists.mean(axis=0)
    mass_2_dist_ci = 2 * mass_2_dists.std(axis=0)

    num_dist_mean = num_dists.mean(axis=0)
    mass_dist_mean = mass_dists.mean(axis=0)

    stats = numpy.zeros(8)
    stats[0] = num_1_err_mean
    stats[1] = num_1_err_ci
    stats[2] = num_2_err_mean
    stats[3] = num_2_err_ci
    stats[4] = mass_1_err_mean
    stats[5] = mass_1_err_ci
    stats[6] = mass_2_err_mean
    stats[7] = mass_2_err_ci

    stats_filename = os.path.join(dirname, "stats.txt")
    numpy.savetxt(stats_filename, stats)

    numpy.savetxt(os.path.join(dirname, "diam.txt"), diam_grid.centers())

    numpy.savetxt(os.path.join(dirname, "num_sect.txt"), aero_binned_single.num_conc)
    numpy.savetxt(os.path.join(dirname, "mass_sect.txt"), aero_binned_single.mass_conc())

    numpy.savetxt(os.path.join(dirname, "num_1_sect.txt"), aero_binned_array.num_conc[0,:])
    numpy.savetxt(os.path.join(dirname, "num_2_sect.txt"), aero_binned_array.num_conc[1,:])
    numpy.savetxt(os.path.join(dirname, "mass_1_sect.txt"), aero_binned_array.mass_conc()[0,:])
    numpy.savetxt(os.path.join(dirname, "mass_2_sect.txt"), aero_binned_array.mass_conc()[1,:])

    numpy.savetxt(os.path.join(dirname, "num_dist_mean.txt"), num_dist_mean)
    numpy.savetxt(os.path.join(dirname, "mass_dist_mean.txt"), mass_dist_mean)

    numpy.savetxt(os.path.join(dirname, "num_1_dist_mean.txt"), num_1_dist_mean)
    numpy.savetxt(os.path.join(dirname, "num_1_dist_ci.txt"), num_1_dist_ci)
    numpy.savetxt(os.path.join(dirname, "num_2_dist_mean.txt"), num_2_dist_mean)
    numpy.savetxt(os.path.join(dirname, "num_2_dist_ci.txt"), num_2_dist_ci)
    numpy.savetxt(os.path.join(dirname, "mass_1_dist_mean.txt"), mass_1_dist_mean)
    numpy.savetxt(os.path.join(dirname, "mass_1_dist_ci.txt"), mass_1_dist_ci)
    numpy.savetxt(os.path.join(dirname, "mass_2_dist_mean.txt"), mass_2_dist_mean)
    numpy.savetxt(os.path.join(dirname, "mass_2_dist_ci.txt"), mass_2_dist_ci)
