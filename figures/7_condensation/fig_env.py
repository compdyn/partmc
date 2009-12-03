#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

maximum_ss = np.zeros([4,4])
 
def make_plot(in_dir, in_file_pattern, out_filename, title, max_ss_i, max_ss_j):
    print in_dir, in_file_pattern
    
    env_state_history = pmc_data_nc.read_history(pmc_data_nc.env_state_t, in_dir, in_file_pattern)
    time = [env_state_history[i][0] for i in range(len(env_state_history))]
    rh = [env_state_history[i][1].relative_humidity for i in range(len(env_state_history))]
    temp = [env_state_history[i][1].temperature for i in range(len(env_state_history))]
    
    print title, (max(rh) - 1)*100.
    
    maximum_ss[max_ss_i, max_ss_j] = (max(rh) - 1)*100.

    plt.clf()
    plt.plot(time,rh,'r')
    ax1 = plt.gca()
    ax1.set_ylim(1,1.003)
    ax2 = plt.twinx()
    ax2.plot(time,temp)
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)


make_plot("../../new_cond/out","cond_1_0001_.*.nc","figs/env_01.pdf","1 hours", 0, 0)
make_plot("../../new_cond/out","cond_2_0001_.*.nc","figs/env_07.pdf","7 hours", 0, 1)
make_plot("../../new_cond/out","cond_3_0001_.*.nc","figs/env_15.pdf","15 hours", 0, 2)
make_plot("../../new_cond/out","cond_4_0001_.*.nc","figs/env_24.pdf","24 hours", 0, 3)

make_plot("../../new_cond/out_comp","cond_1_0001_.*.nc","figs/env_comp_01.pdf","1 hours", 1, 0)
make_plot("../../new_cond/out_comp","cond_2_0001_.*.nc","figs/env_comp_07.pdf","7 hours", 1, 1)
make_plot("../../new_cond/out_comp","cond_3_0001_.*.nc","figs/env_comp_15.pdf","15 hours", 1, 2)
make_plot("../../new_cond/out_comp","cond_4_0001_.*.nc","figs/env_comp_24.pdf","24 hours", 1, 3)

make_plot("../../new_cond/out_size","cond_1_0001_.*.nc","figs/env_size_01.pdf","1 hours", 2, 0)
make_plot("../../new_cond/out_size","cond_2_0001_.*.nc","figs/env_size_07.pdf","7 hours", 2, 1)
make_plot("../../new_cond/out_size","cond_3_0001_.*.nc","figs/env_size_15.pdf","15 hours", 2, 2)
make_plot("../../new_cond/out_size","cond_4_0001_.*.nc","figs/env_size_24.pdf","24 hours", 2, 3)

make_plot("../../new_cond/out_both","cond_1_0001_.*.nc","figs/env_both_01.pdf","1 hours", 3, 0)
make_plot("../../new_cond/out_both","cond_2_0001_.*.nc","figs/env_both_07.pdf","7 hours", 3, 1)
make_plot("../../new_cond/out_both","cond_3_0001_.*.nc","figs/env_both_15.pdf","15 hours", 3, 2)
make_plot("../../new_cond/out_both","cond_4_0001_.*.nc","figs/env_both_24.pdf","24 hours", 3, 3)

np.savetxt("data/maximum_ss.txt", maximum_ss)
