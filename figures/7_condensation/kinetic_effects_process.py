#!/usr/bin/env python

import scipy.io
import sys
import numpy as np

sys.path.append("../../tool")
import partmc
import matplotlib
import mpl_helper
import config
#matplotlib.use("PDF")
#import matplotlib.pyplot as plt


def make_plot(filename_in_080, filename_in_100, filename_in_130, filename_in_150, dir_cloud, in_file_pattern):

## calculate critical supersaturation for each particle
    print filename_in_100

    ncf = scipy.io.netcdf.netcdf_file(filename_in_080, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    particles.sort_by_id()
    env_state_080s = partmc.env_state_t(ncf)
    ncf.close()

    ncf = scipy.io.netcdf.netcdf_file(filename_in_100, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    particles.sort_by_id()
    env_state_100s = partmc.env_state_t(ncf)
    ncf.close()

    dry_diameters = particles.dry_diameters()

    ncf = scipy.io.netcdf.netcdf_file(filename_in_130, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    particles.sort_by_id()
    env_state_130s = partmc.env_state_t(ncf)
    ncf.close()

    ncf = scipy.io.netcdf.netcdf_file(filename_in_150, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    particles.sort_by_id()
    env_state_150s = partmc.env_state_t(ncf)
    ncf.close()

    dry_diameters = particles.dry_diameters()

## calculate time series for RH and maximum RH that occurs
    print dir_cloud, in_file_pattern
    env_state_history = partmc.read_history(partmc.env_state_t, dir_cloud, in_file_pattern)
    env_state_init = env_state_history[0][1]
    time = [env_state_history[i][0] for i in range(len(env_state_history))]
    rh = [env_state_history[i][1].relative_humidity for i in range(len(env_state_history))]
    temp = [env_state_history[i][1].temperature for i in range(len(env_state_history))]

    maximum_ss = (max(rh) - 1) * 100.
    print 'maximum ss ', maximum_ss
    final_ss = (rh[-1] - 1) * 100.  # takes the last value of the time series
    final_rh = rh[-1]
    print 'final ss ', final_ss
 
## calculate time series for each particle
    print dir_cloud, in_file_pattern
    time_filename_list = partmc.get_time_filename_list(dir_cloud, in_file_pattern)
    rh_c = np.zeros((len(dry_diameters),len(time_filename_list)))
    final_rh_c = np.zeros(len(dry_diameters))
    d_c = np.zeros((len(dry_diameters),len(time_filename_list)))
    d = np.zeros((len(dry_diameters),len(time_filename_list)))
    rh_eq = np.zeros((len(dry_diameters),len(time_filename_list)))

    kappas = np.zeros((len(dry_diameters),len(time_filename_list)))
    dry_diam = np.zeros((len(dry_diameters),len(time_filename_list)))
    ids = np.zeros((len(dry_diameters),len(time_filename_list)), dtype = int)
    seconds = np.zeros(len(time_filename_list))
    i_count = 0

    for [time, filename, key] in time_filename_list:
        print time, filename, key
        ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        particles.sort_by_id()
        env_state = partmc.env_state_t(ncf)
        ncf.close()
        rh_c[:,i_count] = particles.critical_rel_humids(env_state)
        d_c[:,i_count] = particles.critical_diameters(env_state)
        d[:,i_count] = particles.diameters()
        rh_eq[:,i_count] = particles.equilib_rel_humids(env_state)
        kappas[:,i_count] = particles.kappas()
        dry_diam[:,i_count] = particles.dry_diameters()

        ids[:,i_count] = particles.ids
        seconds[i_count] = i_count
        i_count = i_count + 1

    d_max = d.max(axis = 1)
    final_rh_c[:] = rh_c[:,-1]

## find the particles in the four categories
## 1. particles that do not activate

    not_activate = np.logical_and(np.all(rh < rh_c, axis=1), np.all(d<d_c, axis=1))

    id_list_not_activate = particles.ids[not_activate]
    print "id_list_not_activate ", id_list_not_activate
    plot_id = id_list_not_activate[0]
    plot_index = partmc.find_nearest_index(particles.ids, plot_id)
    diam_not_activate = d[plot_index,:]
    rh_eq_not_activate = rh_eq[plot_index,:]
    ids_not_activate = ids[plot_index,:]
    rh_c_not_activate = rh_c[plot_index,:]
    d_c_not_activate = d_c[plot_index,:]

    dry_diam_not_activate = dry_diam[plot_index,0]
    kappa_not_activate = kappas[plot_index,0]
    wet_diameters_not_activate = partmc.log_grid(min=1.1*dry_diam_not_activate, max=100*dry_diam_not_activate, n_bin=100).edges()

    equilib_rhs_not_activate_grid_080 = partmc.equilib_rel_humids(env_state_080s,
            kappa_not_activate, dry_diam_not_activate, wet_diameters_not_activate)
    equilib_rhs_not_activate_grid_100 = partmc.equilib_rel_humids(env_state_100s,
            kappa_not_activate, dry_diam_not_activate, wet_diameters_not_activate)
    equilib_rhs_not_activate_grid_130 = partmc.equilib_rel_humids(env_state_130s,
            kappa_not_activate, dry_diam_not_activate, wet_diameters_not_activate)
    equilib_rhs_not_activate_grid_150 = partmc.equilib_rel_humids(env_state_150s,
            kappa_not_activate, dry_diam_not_activate, wet_diameters_not_activate)

## 2. evaporation type     
  
    evaporation = np.logical_and(np.logical_and(np.any(rh > rh_c, axis=1), np.all(d < d_c, axis=1)), (final_rh < final_rh_c))

    id_list_evaporation = particles.ids[evaporation]
    print "id_list_evaporation ", id_list_evaporation

#    plot_id = id_list_evaporation[0]
    rh_c_init = rh_c[:,0]
    rh_c_init_evaporation = np.ma.masked_where(np.logical_not(evaporation), rh_c_init)
    plot_id = rh_c_init_evaporation.argmax()

    plot_index = partmc.find_nearest_index(particles.ids, plot_id)
    diam_evaporation = d[plot_index,:]
    ids_evaporation = ids[plot_index,:]
    rh_eq_evaporation = rh_eq[plot_index,:]
    rh_c_evaporation = rh_c[plot_index,:]
    d_c_evaporation = d_c[plot_index,:]

    dry_diam_evaporation = dry_diam[plot_index,0]
    kappa_evaporation = kappas[plot_index,0]
    wet_diameters_evaporation = partmc.log_grid(min=1.1*dry_diam_evaporation, max=100*dry_diam_evaporation, n_bin=100).edges()

    equilib_rhs_evaporation_grid_080 = partmc.equilib_rel_humids(env_state_080s,
            kappa_evaporation, dry_diam_evaporation, wet_diameters_evaporation)
    equilib_rhs_evaporation_grid_100 = partmc.equilib_rel_humids(env_state_100s,
            kappa_evaporation, dry_diam_evaporation, wet_diameters_evaporation)
    equilib_rhs_evaporation_grid_130 = partmc.equilib_rel_humids(env_state_130s,
            kappa_evaporation, dry_diam_evaporation, wet_diameters_evaporation)
    equilib_rhs_evaporation_grid_150 = partmc.equilib_rel_humids(env_state_150s,
            kappa_evaporation, dry_diam_evaporation, wet_diameters_evaporation)

## 3. deactivation type

    deactivation = np.logical_and(np.any(rh > rh_c, axis=1), np.any(d > d_c, axis=1))

    id_list_deactivation = particles.ids[deactivation]
    print "id_list_deactivation ", id_list_deactivation

#    plot_id = id_list_deactivation[0]

    rh_c_init = rh_c[:,0]
    rh_c_init_deactivation = np.ma.masked_where(np.logical_not(deactivation), rh_c_init)
    plot_id = rh_c_init_deactivation.argmax()

    plot_index = partmc.find_nearest_index(particles.ids, plot_id)
    diam_deactivation = d[plot_index,:]
    ids_deactivation = ids[plot_index,:]
    rh_eq_deactivation = rh_eq[plot_index,:]
    rh_c_deactivation = rh_c[plot_index,:]
    d_c_deactivation = d_c[plot_index,:]

    dry_diam_deactivation = dry_diam[plot_index,0]
    kappa_deactivation = kappas[plot_index,0]
    wet_diameters_deactivation = partmc.log_grid(min=1.1*dry_diam_deactivation, max=100*dry_diam_deactivation, n_bin=100).edges()

    equilib_rhs_deactivation_grid_080 = partmc.equilib_rel_humids(env_state_080s,
            kappa_deactivation, dry_diam_deactivation, wet_diameters_deactivation)
    equilib_rhs_deactivation_grid_100 = partmc.equilib_rel_humids(env_state_100s,
            kappa_deactivation, dry_diam_deactivation, wet_diameters_deactivation)
    equilib_rhs_deactivation_grid_130 = partmc.equilib_rel_humids(env_state_130s,
            kappa_deactivation, dry_diam_deactivation, wet_diameters_deactivation)
    equilib_rhs_deactivation_grid_150 = partmc.equilib_rel_humids(env_state_150s,
            kappa_deactivation, dry_diam_deactivation, wet_diameters_deactivation)

## 4. inertial type

    inertial = np.logical_and(np.logical_and(np.any(rh > rh_c, axis=1), np.all(d < d_c, axis=1)), (final_rh > final_rh_c))

    id_list_inertial = particles.ids[inertial] 
    print "id_list_inertial ", id_list_inertial

    plot_id = id_list_inertial[0]
    plot_index = partmc.find_nearest_index(particles.ids, plot_id)
    diam_inertial = d[plot_index,:]
    ids_inertial = ids[plot_index,:]
    rh_eq_inertial = rh_eq[plot_index,:]
    rh_c_inertial = rh_c[plot_index,:]
    d_c_inertial = d_c[plot_index,:]

    dry_diam_inertial = dry_diam[plot_index,0]
    kappa_inertial = kappas[plot_index,0]
    wet_diameters_inertial = partmc.log_grid(min=1.1*dry_diam_inertial, max=100*dry_diam_inertial, n_bin=100).edges()

    equilib_rhs_inertial_grid_080 = partmc.equilib_rel_humids(env_state_080s,
            kappa_inertial, dry_diam_inertial, wet_diameters_inertial)
    equilib_rhs_inertial_grid_100 = partmc.equilib_rel_humids(env_state_100s,
            kappa_inertial, dry_diam_inertial, wet_diameters_inertial)
    equilib_rhs_inertial_grid_130 = partmc.equilib_rel_humids(env_state_130s,
            kappa_inertial, dry_diam_inertial, wet_diameters_inertial)
    equilib_rhs_inertial_grid_150 = partmc.equilib_rel_humids(env_state_150s,
            kappa_inertial, dry_diam_inertial, wet_diameters_inertial)

    np.savetxt("seconds.txt", seconds)
    np.savetxt("rh.txt", rh)

    np.savetxt("diam_not_activate.txt", diam_not_activate)
    np.savetxt("diam_evaporation.txt", diam_evaporation)
    np.savetxt("diam_deactivation.txt", diam_deactivation)
    np.savetxt("diam_inertial.txt", diam_inertial)

    np.savetxt("rh_eq_not_activate.txt", rh_eq_not_activate)
    np.savetxt("rh_eq_evaporation.txt", rh_eq_evaporation)
    np.savetxt("rh_eq_deactivation.txt", rh_eq_deactivation)
    np.savetxt("rh_eq_inertial.txt", rh_eq_inertial)

    np.savetxt("ids_not_activate.txt", ids_not_activate)
    np.savetxt("ids_evaporation.txt", ids_evaporation)
    np.savetxt("ids_deactivation.txt", ids_deactivation)
    np.savetxt("ids_inertial.txt", ids_inertial)

    np.savetxt("rh_c_not_activate.txt", rh_c_not_activate)
    np.savetxt("rh_c_evaporation.txt", rh_c_evaporation)
    np.savetxt("rh_c_deactivation.txt", rh_c_deactivation)
    np.savetxt("rh_c_inertial.txt", rh_c_inertial)

    np.savetxt("d_c_not_activate.txt", d_c_not_activate)
    np.savetxt("d_c_evaporation.txt", d_c_evaporation)
    np.savetxt("d_c_deactivation.txt", d_c_deactivation)
    np.savetxt("d_c_inertial.txt", d_c_inertial)

    np.savetxt("wet_diameters_not_activate.txt", wet_diameters_not_activate)
    np.savetxt("wet_diameters_evaporation.txt", wet_diameters_evaporation)
    np.savetxt("wet_diameters_deactivation.txt", wet_diameters_deactivation)
    np.savetxt("wet_diameters_inertial.txt", wet_diameters_inertial)

    np.savetxt("equilib_rhs_not_activate_grid_080.txt", equilib_rhs_not_activate_grid_080)
    np.savetxt("equilib_rhs_evaporation_grid_080.txt", equilib_rhs_evaporation_grid_080)
    np.savetxt("equilib_rhs_deactivation_grid_080.txt", equilib_rhs_deactivation_grid_080)
    np.savetxt("equilib_rhs_inertial_grid_080.txt", equilib_rhs_inertial_grid_080)

    np.savetxt("equilib_rhs_not_activate_grid_100.txt", equilib_rhs_not_activate_grid_100)
    np.savetxt("equilib_rhs_evaporation_grid_100.txt", equilib_rhs_evaporation_grid_100)
    np.savetxt("equilib_rhs_deactivation_grid_100.txt", equilib_rhs_deactivation_grid_100)
    np.savetxt("equilib_rhs_inertial_grid_100.txt", equilib_rhs_inertial_grid_100)

    np.savetxt("equilib_rhs_not_activate_grid_130.txt", equilib_rhs_not_activate_grid_130)
    np.savetxt("equilib_rhs_evaporation_grid_130.txt", equilib_rhs_evaporation_grid_130)
    np.savetxt("equilib_rhs_deactivation_grid_130.txt", equilib_rhs_deactivation_grid_130)
    np.savetxt("equilib_rhs_inertial_grid_130.txt", equilib_rhs_inertial_grid_130)

    np.savetxt("equilib_rhs_not_activate_grid_150.txt", equilib_rhs_not_activate_grid_150)
    np.savetxt("equilib_rhs_evaporation_grid_150.txt", equilib_rhs_evaporation_grid_150)
    np.savetxt("equilib_rhs_deactivation_grid_150.txt", equilib_rhs_deactivation_grid_150)
    np.savetxt("equilib_rhs_inertial_grid_150.txt", equilib_rhs_inertial_grid_150)

filename_in_080 = "condensation_data/cond_tenthou_13_ref_ver210_mem_10_0001_00000081.nc"
filename_in_100 = "condensation_data/cond_tenthou_13_ref_ver210_mem_10_0001_00000101.nc"
filename_in_130 = "condensation_data/cond_tenthou_13_ref_ver210_mem_10_0001_00000131.nc"
filename_in_150 = "condensation_data/cond_tenthou_13_ref_ver210_mem_10_0001_00000151.nc"

dir_cloud = "condensation_data"
in_file_pattern = "cond_tenthou_13_ref_ver210_mem_10_0001_.*.nc"

make_plot(filename_in_080, filename_in_100, filename_in_130, filename_in_150, dir_cloud, in_file_pattern)
