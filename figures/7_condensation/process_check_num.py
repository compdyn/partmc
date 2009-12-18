#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
const = pmc_data_nc.load_constants("../../src/constants.f90")

def check_num(in_dir1, in_filename1, in_file_pattern, indir2, in_filename2, out_filename1, out_filename2, out_filename3, counter):
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir1+in_filename1)
    particles1 = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir2+in_filename2)
    particles2 = pmc_data_nc.aero_particle_array_t(ncf)
    particles2.aero_data.kappa[17] = 0.1
    particles2.aero_data = particles1.aero_data
    ncf.close()

    final_wet_diameter = particles1.diameter()
    is_activated1 = (final_wet_diameter > 3e-6)
    sum_tot1 = sum(1/particles1.comp_vol) * particles1.comp_vol[0]
    num_act1 = sum(1/particles1.comp_vol[is_activated1]) * particles1.comp_vol[0] 
    id_list_act1 = particles1.id[is_activated1]
    is_not_activated1 = np.logical_not(is_activated1)
    id_list_not_act1 = particles1.id[is_not_activated1] 
    ccn_cn_ratio1 = num_act1 / sum_tot1

    env_state_history = pmc_data_nc.read_history(pmc_data_nc.env_state_t, in_dir1, in_file_pattern)
    time = [env_state_history[i][0] for i in range(len(env_state_history))]
    rh = [env_state_history[i][1].relative_humidity for i in range(len(env_state_history))]
    maximum_ss = (max(rh) - 1)*100.
    max_index = np.argmax(np.array(rh))
    time_index = time[max_index]    
    time_filename_list = pmc_data_nc.get_time_filename_list(in_dir1, in_file_pattern)
    max_filename = pmc_data_nc.find_filename_at_time(time_filename_list, time_index)
    ncf = Scientific.IO.NetCDF.NetCDFFile(max_filename)
    max_env_state = pmc_data_nc.env_state_t(ncf)
    ncf.close()


    max_wet_diameter = np.zeros_like(final_wet_diameter)
    max_critical_ratio = np.zeros_like(final_wet_diameter)
    
    for [time, filename, key] in time_filename_list:
        print 'time filename key ', time, filename, key
        ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
        particles = pmc_data_nc.aero_particle_array_t(ncf) 
        env_state = pmc_data_nc.env_state_t(ncf)
        ncf.close()
        
        wet_diameter = particles.diameter()
        max_wet_diameter = np.maximum(max_wet_diameter, wet_diameter)
        critical_diameter = particles.critical_diameter(env_state, const)
        critical_ratio = wet_diameter / critical_diameter
        max_critical_ratio = np.maximum(max_critical_ratio, critical_ratio)

    max_wet_ratio = max_wet_diameter / final_wet_diameter
    
    s_crit = (particles2.critical_rh(max_env_state, const) - 1)*100
    is_activated2 = (s_crit <= maximum_ss)
    id_list_act2 = particles2.id[is_activated2]
    is_not_activated2 = np.logical_not(is_activated2)
    id_list_not_act2 = particles2.id[is_not_activated2]

    sum_tot2 = sum(1/particles2.comp_vol) * particles2.comp_vol[0]
    num_act2 = sum(1/particles2.comp_vol[is_activated2]) * particles2.comp_vol[0]
    ccn_cn_ratio2 = num_act2 / sum_tot2

    set_act1 = set(id_list_act1)
    set_not_act1 = set(id_list_not_act1)
    set_act2 = set(id_list_act2)
    set_not_act2 = set(id_list_not_act2)

    act_1_2 = (set_act1 & set_act2)
    not_1_2 = (set_not_act1 &  set_not_act2)
    act1_not_act2 = (set_act1 & set_not_act2)
    act2_not_act1 = (set_act2 & set_not_act1)

    print 'act_1_2 ', act_1_2
    print 'not_1_2 ', not_1_2
    print 'act1_not_act2 ', act1_not_act2
    print 'act2_not_act1 ', act2_not_act1

    diam_by_id1 = dict(zip(particles1.id, particles1.dry_diameter()))
    diam_by_id2 = dict(zip(particles2.id, particles2.dry_diameter()))
    scrit_by_id1 = dict(zip(particles1.id, (particles1.critical_rh(max_env_state,const) - 1)*100))
    scrit_by_id2 = dict(zip(particles2.id, (particles2.critical_rh(max_env_state,const) - 1)*100))
    oc_by_id1 = dict(zip(particles1.id, particles1.mass(include = ["BC"])/particles1.mass(exclude=["H2O"])))
    oc_by_id2 = dict(zip(particles2.id, particles2.mass(include = ["BC"])/particles2.mass(exclude=["H2O"])))
    wet_ratio_by_id1 = dict(zip(particles1.id, max_wet_ratio)) 
    critical_ratio_by_id1 = dict(zip(particles1.id, max_critical_ratio))

    diam_act_1_2 = [diam_by_id2[id] for id in act_1_2]
    scrit_act_1_2 = [scrit_by_id2[id] for id in act_1_2]
    oc_act_1_2 = [oc_by_id2[id] for id in act_1_2]
    wet_ratio_act_1_2 = [wet_ratio_by_id1[id] for id in act_1_2]
    critical_ratio_act_1_2 = [critical_ratio_by_id1[id] for id in act_1_2]

    diam_not_1_2 = [diam_by_id2[id] for id in not_1_2]
    scrit_not_1_2 = [scrit_by_id2[id] for id in not_1_2]
    oc_not_1_2 = [oc_by_id2[id] for id in not_1_2]
    wet_ratio_not_1_2 = [wet_ratio_by_id1[id] for id in not_1_2]
    critical_ratio_not_1_2 = [critical_ratio_by_id1[id] for id in not_1_2]
    
    diam_act1_not_act2 = [diam_by_id2[id] for id in act1_not_act2]
    scrit_act1_not_act2 = [scrit_by_id2[id] for id in act1_not_act2]
    oc_act1_not_act2 = [oc_by_id2[id] for id in act1_not_act2]
    wet_ratio_act1_not_act2 = [wet_ratio_by_id1[id] for id in act1_not_act2]
    critical_ratio_act1_not_act2 = [critical_ratio_by_id1[id] for id in act1_not_act2]
    
    diam_act2_not_act1 = [diam_by_id2[id] for id in act2_not_act1]
    scrit_act2_not_act1 = [scrit_by_id2[id] for id in act2_not_act1]
    oc_act2_not_act1 = [oc_by_id2[id] for id in act2_not_act1]
    wet_ratio_act2_not_act1 = [wet_ratio_by_id1[id] for id in act2_not_act1]
    critical_ratio_act2_not_act1 = [critical_ratio_by_id1[id] for id in act2_not_act1]

    plt.figure()
    plt.semilogy(seconds, d[0,:], 'b-', label = 'act1_2')
    plt.semilogy(seconds, d[1,:], 'g-', label = 'not_1_2')
    plt.semilogy(seconds, d[2,:], 'r-', label = 'act2_not_act1')
    plt.semilogy(seconds, d[3,:], 'b-', label = 'act1_2')
    plt.semilogy(seconds, d[4,:], 'g-', label = 'not_1_2')
    plt.semilogy(seconds, d[5,:], 'r-', label = 'act2_not_act1')
    plt.xlabel("time (s)")
    plt.ylabel("diameter (m)")
    plt.legend(loc = 'upper left')
    fig = plt.gcf()
    fig.savefig('diameters.pdf')
    
    plt.figure()
    plt.loglog(diam_act_1_2, scrit_act_1_2, 'bx')
    plt.loglog(diam_not_1_2, scrit_not_1_2, 'gx') 
    plt.loglog(diam_act1_not_act2, scrit_act1_not_act2, 'kx')
    plt.loglog(diam_act2_not_act1, scrit_act2_not_act1, 'rx')
    plt.axhline(y=maximum_ss)
    plt.xlabel("diameter (m)")
    plt.ylabel("critical supersaturation (%)")
    fig = plt.gcf()
    fig.savefig(out_filename1)

    plt.figure()
    plt.semilogx(diam_act_1_2, wet_ratio_act_1_2, 'bx')
    plt.semilogx(diam_not_1_2, wet_ratio_not_1_2, 'gx')
    plt.semilogx(diam_act1_not_act2, wet_ratio_act1_not_act2, 'kx')
    plt.semilogx(diam_act2_not_act1, wet_ratio_act2_not_act1, 'rx')
    plt.xlabel("diameter (m)")
    plt.ylabel("wet ratio")
    fig = plt.gcf()
    fig.savefig(out_filename2)

    plt.figure()
    plt.loglog(diam_act_1_2, critical_ratio_act_1_2, 'bx')
    plt.loglog(diam_not_1_2, critical_ratio_not_1_2, 'gx')
    plt.loglog(diam_act1_not_act2, critical_ratio_act1_not_act2, 'kx')
    plt.loglog(diam_act2_not_act1, critical_ratio_act2_not_act1, 'rx')
    plt.xlabel("diameter (m)")
    plt.ylabel("critial ratio")
    fig = plt.gcf()
    fig.savefig(out_filename3)

for counter in range(5,6):

    in_dir1 = "../../new_cond/out/"
    in_dir2 = "../../new_cond/start/"
    filename_in1 = "cond_%02d_ref_0001_00000601.nc" % counter
    filename_in2 = "urban_plume_wc_0001_000000%02d.nc" % counter
    in_file_pattern = "cond_%02d_ref_0001_(.*).nc" % counter
    out_filename1 = "figs/check_num_scrit_%02d.pdf" % counter
    out_filename2 = "figs/check_num_wr_%02d.pdf" % counter
    out_filename3 = "figs/check_num_cr_%02d.pdf" % counter
    check_num(in_dir1, filename_in1, in_file_pattern, in_dir2, filename_in2, out_filename1, out_filename2, out_filename3, counter)

