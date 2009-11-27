#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
import config
const = pmc_data_nc.load_constants("../../src/constants.f90")

def make_plot(netcdf_dir, netcdf_pattern, out_filename):
    time_filename_list = pmc_data_nc.get_time_filename_list(netcdf_dir, netcdf_pattern)
    ccn_array = np.zeros([len(time_filename_list),4])
    i_counter = 0
    for [time, filename, key] in time_filename_list:
        print time, filename, key
        ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
        particles = pmc_data_nc.aero_particle_array_t(ncf)
        env_state = pmc_data_nc.env_state_t(ncf)
        ncf.close()

        s_crit = (particles.kappa_rh(env_state, const) - 1)*100

        activated_1 = (s_crit < config.s_crit_1)
        number_act_1 = sum(1/particles.comp_vol[activated_1])

        activated_2 = (s_crit < config.s_crit_2)
        number_act_2 = sum(1/particles.comp_vol[activated_2])

        activated_3 = (s_crit < config.s_crit_3)
        number_act_3 = sum(1/particles.comp_vol[activated_3])

        ccn_array[i_counter,0]= time
        ccn_array[i_counter,1]= number_act_1
        ccn_array[i_counter,2]= number_act_2
        ccn_array[i_counter,3]= number_act_3
        i_counter += 1

    print ccn_array

    np.savetxt(out_filename, ccn_array)

make_plot("../../urban_plume2/out/", "urban_plume_wc_0001_(.*).nc", "data/ccn_conc_wc.txt")
make_plot("../../urban_plume2/out/", "urban_plume_nc_0001_(.*).nc", "data/ccn_conc_nc.txt")



