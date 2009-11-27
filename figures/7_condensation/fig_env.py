#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

def make_plot(in_dir, in_file_pattern, out_filename,title):
    print in_dir, in_file_pattern
    
    env_state_history = pmc_data_nc.read_history(pmc_data_nc.env_state_t, in_dir, in_file_pattern)
    time = [env_state_history[i][0] for i in range(len(env_state_history))]
    rh = [env_state_history[i][1].relative_humidity for i in range(len(env_state_history))]
    temp = [env_state_history[i][1].temperature for i in range(len(env_state_history))]
    
    print title, (max(rh) - 1)*100.

    plt.clf()
    plt.plot(time,rh,'r')
    ax1 = plt.gca()
    ax1.set_ylim(1,1.003)
    ax2 = plt.twinx()
    ax2.plot(time,temp)
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

make_plot("../../new_cond/out","cond_1_0001_.*.nc","figs/env_01.pdf","1 hours")
make_plot("../../new_cond/out","cond_2_0001_.*.nc","figs/env_07.pdf","7 hours")
make_plot("../../new_cond/out","cond_3_0001_.*.nc","figs/env_15.pdf","15 hours")
make_plot("../../new_cond/out","cond_4_0001_.*.nc","figs/env_24.pdf","24 hours")

