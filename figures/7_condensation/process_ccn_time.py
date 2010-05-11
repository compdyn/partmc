#!/usr/bin/env python2.5

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import config

netcdf_dir = "../../scenarios/2_urban_plume2/out/"
netcdf_pattern = "urban_plume_wc_0001_(.*).nc"
time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

ccn_cn_array = np.zeros([len(time_filename_list),4])
i_counter = 0
for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    total_number = len(particles.ids)
    s_crit = (particles.critical_rel_humids(env_state) - 1)*100
    
    activated_1 = (s_crit < config.s_crit_1)
    number_act_1 = sum(activated_1)
    ccn_cn_1 = float(number_act_1) / total_number

    activated_2 = (s_crit < config.s_crit_2)
    number_act_2 = sum(activated_2)
    ccn_cn_2 = float(number_act_2) / total_number

    activated_3 = (s_crit < config.s_crit_3)
    number_act_3 = sum(activated_3)
    ccn_cn_3 = float(number_act_3) / total_number

    ccn_cn_array[i_counter,0]= time
    ccn_cn_array[i_counter,1]= ccn_cn_1
    ccn_cn_array[i_counter,2]= ccn_cn_2
    ccn_cn_array[i_counter,3]= ccn_cn_3
    i_counter += 1

print ccn_cn_array

np.savetxt("data/ccn_cn.txt", ccn_cn_array)

    



