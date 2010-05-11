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

def make_plot(netcdf_dir, netcdf_pattern, out_filename):
    time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)
    ccn_array = np.zeros([len(time_filename_list),4])
    i_counter = 0
    for [time, filename, key] in time_filename_list:
        print time, filename, key
        ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        s_crit = (particles.critical_rel_humids(env_state) - 1)*100

        activated_1 = (s_crit < config.s_crit_1)
        number_act_1 = sum(1/particles.comp_vols[activated_1])

        activated_2 = (s_crit < config.s_crit_2)
        number_act_2 = sum(1/particles.comp_vols[activated_2])

        activated_3 = (s_crit < config.s_crit_3)
        number_act_3 = sum(1/particles.comp_vols[activated_3])

        ccn_array[i_counter,0]= time
        ccn_array[i_counter,1]= number_act_1
        ccn_array[i_counter,2]= number_act_2
        ccn_array[i_counter,3]= number_act_3
        i_counter += 1

    print ccn_array

    np.savetxt(out_filename, ccn_array)

make_plot("../../scenarios/2_urban_plume2/out/", "urban_plume_wc_0001_(.*).nc", "data/ccn_conc_wc.txt")
make_plot("../../scenarios/2_urban_plume2/out/", "urban_plume_nc_0001_(.*).nc", "data/ccn_conc_nc.txt")



