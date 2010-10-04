#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import config

def make_plot(hour_counter, case_counter, in_files):
    x_axis = partmc.log_grid(min=1e-9,max=1e-5,n_bin=70)
    y_axis = partmc.log_grid(min=1e-3,max=1e2,n_bin=50)
    x_centers = x_axis.centers()
    y_centers = y_axis.centers()
    i_counter = 0
    hist_array = np.zeros([len(x_centers), len(y_centers), config.i_loop_max])
    hist_average = np.zeros([len(x_centers), len(y_centers)])
    hist_std = np.zeros([len(x_centers), len(y_centers)])
    hist_std_norm = np.zeros([len(x_centers), len(y_centers)])
    ccn_array = np.zeros([4,config.i_loop_max])  # s_crit values x i_loop
    bc_array = np.zeros([4,config.i_loop_max])

    for file in in_files:
        ncf = scipy.io.netcdf.netcdf_file(config.netcdf_dir+'/'+file, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        dry_diameters = particles.dry_diameters()
        total_number = sum(1/particles.comp_vols)
        bc = particles.masses(include = ["BC"])
        total_bc = sum(bc/particles.comp_vols)

        s_crit = (particles.critical_rel_humids(env_state) - 1)*100

        activated_1 = (s_crit < config.s_crit_1)
        number_act_1 = sum(1/particles.comp_vols[activated_1])
        bc_mass_act1 = sum(bc[activated_1]/particles.comp_vols[activated_1])

        activated_2 = (s_crit < config.s_crit_2)
        number_act_2 = sum(1/particles.comp_vols[activated_2])
        bc_mass_act2 = sum(bc[activated_2]/particles.comp_vols[activated_2])

        activated_3 = (s_crit < config.s_crit_3)
        number_act_3 = sum(1/particles.comp_vols[activated_3])
        bc_mass_act3 = sum(bc[activated_3]/particles.comp_vols[activated_3])

        activated_4 = (s_crit < config.s_crit_4)
        number_act_4 = sum(1/particles.comp_vols[activated_4])
        bc_mass_act4 = sum(bc[activated_4]/particles.comp_vols[activated_4])

        ccn_array[0,i_counter]= float(number_act_1)/total_number
        ccn_array[1,i_counter]= float(number_act_2)/total_number
        ccn_array[2,i_counter]= float(number_act_3)/total_number
        ccn_array[3,i_counter]= float(number_act_4)/total_number

        bc_array[0,i_counter]= float(bc_mass_act1)/total_bc
        bc_array[1,i_counter]= float(bc_mass_act2)/total_bc
        bc_array[2,i_counter]= float(bc_mass_act3)/total_bc
        bc_array[3,i_counter]= float(bc_mass_act4)/total_bc

        i_counter += 1

    print 'counters 3 ', hour_counter, case_counter

    ccn = np.average(ccn_array, axis = 1)
    ccn_average[hour_counter,case_counter,:] = ccn
    ccn_std[hour_counter,case_counter,:] = np.std(ccn_array, axis = 1)

    bc_average[hour_counter,case_counter,:] = np.average(bc_array, axis = 1)
    bc_std[hour_counter, case_counter, :] = np.std(bc_array, axis = 1)

hour_counter = 0
ccn_average = np.zeros([25,21,4])
ccn_std = np.zeros([25,21,4])
ccn_average_overall = np.zeros([21,4])
ccn_std_overall = np.zeros([21,4])

bc_average = np.zeros([25,21,4])
bc_std = np.zeros([25,21,4])
bc_average_overall = np.zeros([21,4])
bc_std_overall = np.zeros([21,4])

for hour in range(1,26):
    print "hour = ", hour
    case_counter = 0
    print 'counters 1 ', hour_counter, case_counter
    for counter in ["1K_wei+1", "1K_flat", "1K_wei-1", "1K_wei-2", "1K_wei-3", "1K_wei-4", 
	 		"10K_wei+1", "10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4",
			"100K_wei+1", "100K_flat", "100K_wei-1", "100K_wei-2", "100K_wei-3", "100K_wei-4",
			"1K_mfa", "10K_mfa", "100K_mfa"]:
        print 'counter ', counter
        files = []
        for i_loop in range(0,config.i_loop_max):
            filename_in = "urban_plume_wc_%s_0%03d_000000%02d.nc" % (counter,i_loop+1,hour)
            files.append(filename_in)
        make_plot(hour_counter, case_counter, files)
        print 'counters 2', hour_counter, case_counter
        case_counter = case_counter + 1
    hour_counter = hour_counter + 1

ccn_std = ccn_std / ccn_average
bc_std = bc_std / bc_average

ccn_std_overall = np.average(ccn_std, axis = 0)
bc_std_overall = np.average(bc_std, axis = 0)

ccn_average_overall = np.average(ccn_average, axis = 0)
bc_average_overall = np.average(bc_average, axis = 0)

np.savetxt("data/ccn_average_ss1.txt", ccn_average[:,:,0])
np.savetxt("data/ccn_average_ss2.txt", ccn_average[:,:,1])
np.savetxt("data/ccn_average_ss3.txt", ccn_average[:,:,2])
np.savetxt("data/ccn_average_ss4.txt", ccn_average[:,:,3])

np.savetxt("data/bc_average_ss1.txt", bc_average[:,:,0])
np.savetxt("data/bc_average_ss2.txt", bc_average[:,:,1])
np.savetxt("data/bc_average_ss3.txt", bc_average[:,:,2])
np.savetxt("data/bc_average_ss4.txt", bc_average[:,:,3])

np.savetxt("data/ccn_std_ss1.txt", ccn_std[:,:,0])
np.savetxt("data/ccn_std_ss2.txt", ccn_std[:,:,1])
np.savetxt("data/ccn_std_ss3.txt", ccn_std[:,:,2])
np.savetxt("data/ccn_std_ss4.txt", ccn_std[:,:,3])

np.savetxt("data/bc_std_ss1.txt", bc_std[:,:,0])
np.savetxt("data/bc_std_ss2.txt", bc_std[:,:,1])
np.savetxt("data/bc_std_ss3.txt", bc_std[:,:,2])
np.savetxt("data/bc_std_ss4.txt", bc_std[:,:,3])

np.savetxt("data/ccn_std_overall_ss1.txt", ccn_std_overall[:,0])
np.savetxt("data/ccn_std_overall_ss2.txt", ccn_std_overall[:,1])
np.savetxt("data/ccn_std_overall_ss3.txt", ccn_std_overall[:,2])
np.savetxt("data/ccn_std_overall_ss4.txt", ccn_std_overall[:,3])

np.savetxt("data/bc_std_overall_ss1.txt", bc_std_overall[:,0])
np.savetxt("data/bc_std_overall_ss2.txt", bc_std_overall[:,1])
np.savetxt("data/bc_std_overall_ss3.txt", bc_std_overall[:,2])
np.savetxt("data/bc_std_overall_ss4.txt", bc_std_overall[:,3])

np.savetxt("data/ccn_average_overall_ss1.txt", ccn_average_overall[:,0])
np.savetxt("data/ccn_average_overall_ss2.txt", ccn_average_overall[:,1])
np.savetxt("data/ccn_average_overall_ss3.txt", ccn_average_overall[:,2])
np.savetxt("data/ccn_average_overall_ss4.txt", ccn_average_overall[:,3])

np.savetxt("data/bc_average_overall_ss1.txt", bc_average_overall[:,0])
np.savetxt("data/bc_average_overall_ss2.txt", bc_average_overall[:,1])
np.savetxt("data/bc_average_overall_ss3.txt", bc_average_overall[:,2])
np.savetxt("data/bc_average_overall_ss4.txt", bc_average_overall[:,3])

