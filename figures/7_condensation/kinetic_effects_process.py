#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(filename_in, dir_cloud, in_file_pattern):

## calculate critical supersaturation for each particle
    print filename_in
    ncf = scipy.io.netcdf.netcdf_file(filename_in, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state_plume = partmc.env_state_t(ncf)
    ncf.close()

    dry_diameters = particles.dry_diameters()
    s_crit = (particles.critical_rel_humids(env_state_plume) - 1)*100
    d_crit = particles.critical_diameters(env_state_plume)
    kappas = particles.kappas()

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
    print 'final ss ', final_ss
 
## calculate time series of diameters of each particle
    print dir_cloud, in_file_pattern
    time_filename_list = partmc.get_time_filename_list(dir_cloud, in_file_pattern)
    d = np.zeros((len(dry_diameters),len(time_filename_list)))
    seconds = np.zeros(len(time_filename_list))
    i_count = 0

    for [time, filename, key] in time_filename_list:
        print time, filename, key
        ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        ncf.close()
        wet_diameters = particles.diameters()
        d[:,i_count] = wet_diameters
        seconds[i_count] = i_count
        i_count = i_count + 1

    print "diameter for particle 1 ", d[1,:]
 
    d_max = d.max(axis = 1)

    print "dmax ", d_max   
    
## find the particles in the four categories
## 1. particles that do now activate

    not_activate = ((maximum_ss < s_crit) & (d_max < d_crit))
    id_list_not_activate = particles.ids[not_activate]
    print "id_list_not_activate ", id_list_not_activate
    plot_id = id_list_not_activate[0]
    plot_index = partmc.find_nearest_index(particles.ids, plot_id)
    diam_not_activate = d[plot_index,:]
  
    dry_diameter = dry_diameters[plot_index]
    kappa = kappas[plot_index]
    wet_diameters = partmc.log_grid(min = diam_not_activate.min(), max = diam_not_activate.max(), n_bin=100).edges()
    equilib_rhs = partmc.equilib_rel_humids(env_state_init,
            kappa, dry_diameter, wet_diameters)
     
    print "diam_not_activate, rh ", diam_not_activate, rh
    plt.clf()
    plt.plot(diam_not_activate*1e6, rh)
    plt.hold(True)
    plt.plot(wet_diameters*1e6, equilib_rhs, 'r') 
    fig = plt.gcf()
    fig.savefig("not_activate.pdf")


## 2. evaporation type     
  
    evaporation = ((maximum_ss > s_crit) & (d_max < d_crit))
    id_list_evaporation = particles.ids[evaporation]
    print "id_list_evaporation ", id_list_evaporation
## 3. deactivation type

    deactivation = ((maximum_ss > s_crit) & (d_max > d_crit))
    id_list_deactivation = particles.ids[deactivation]
    print "id_list_deactivation ", id_list_deactivation
## 4. inertal type

    inertial = ((maximum_ss > s_crit) & (d_max < d_crit))
    id_list_inertial = particles.ids[inertial] 
    print "id_list_inertial ", id_list_inertial

filename_in = "/home/ching1/subversion/partmc/branches/joseph2/scenarios/3_condense/out/ver210_ref_13/cond_tenthou_13_ref_ver210_mem_10_0001_00000013.nc"
dir_cloud = "/home/ching1/subversion/partmc/branches/joseph2/scenarios/3_condense/out/ver210_ref_13"
in_file_pattern = "cond_tenthou_13_ref_ver210_mem_10_0001_.*.nc"

make_plot(filename_in, dir_cloud, in_file_pattern)
