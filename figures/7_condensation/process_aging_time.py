#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import config

netcdf_dir = "../../scenarios/2_urban_plume2/out/"
netcdf_pattern = "urban_plume_nc_0001_(.*).nc"
time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

ccn_cn_array = np.zeros([len(time_filename_list),4])
i_counter = 0
age_by_id = {}
for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    s_crit = (particles.critical_rel_humids(env_state) - 1)*100
    
    activated_1 = (s_crit < config.s_crit_1)
    activated_id_set = set(particles.ids[activated_1])

    diam_by_id = dict(zip(particles.ids, particles.dry_diameters()))
    time_by_id = dict(zip(particles.ids, particles.least_create_times))
    if (i_counter > 0):
        aged_id = activated_id_set - activated_previous
        for id in aged_id:
            age_by_id[id] = env_state.elapsed_time - time_by_id[id]
    activated_previous = activated_id_set
    i_counter += 1



    



