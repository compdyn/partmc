#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import mpl_helper
import partmc

def make_plot(in_dir1, in_dir2, in_file_pattern, out_filename):
#def make_plot(in_dir1, in_file_pattern, out_filename):
    print in_dir1, in_file_pattern
    

    (figure, axes) = mpl_helper.make_fig(right_margin=0.8)
    env_state_history = partmc.read_history(partmc.env_state_t, in_dir1, in_file_pattern)
    time = [env_state_history[i][0] for i in range(len(env_state_history))]
    rh1 = [env_state_history[i][1].relative_humidity for i in range(len(env_state_history))]
    temp = [env_state_history[i][1].temperature for i in range(len(env_state_history))]
    print time
    
    env_state_history = partmc.read_history(partmc.env_state_t, in_dir2, in_file_pattern)
    time = [env_state_history[i][0] for i in range(len(env_state_history))]
    rh2 = [env_state_history[i][1].relative_humidity for i in range(len(env_state_history))]
    temp = [env_state_history[i][1].temperature for i in range(len(env_state_history))]


    for i in range(len(env_state_history)):
        time[i] = (time[i] - 14 * 3600.)/60

#    print (max(rh1) - 1)*100., (min(rh1) - 1)*100. 
    
    axes.plot(time,rh1,'r', label = 'w. e.')
    axes.plot(time,rh2,'g', label = 'n. e.')

    axes2 = axes.twinx()
    axes.set_ylim(0.7,1.1)
    axes.set_ylabel('RH')
    axes.legend(loc = 'lower center')
    axes2.plot(time,temp)
    axes2.set_ylabel("Temperature / K")
    axes.set_xlabel("time / min")
    axes.grid(True)
    figure.savefig(out_filename)

dir_name1 = "../../scenarios/8_condense_p/out_a"
dir_name2 = "../../scenarios/8_condense_p/out_a_no"

filename_in1 = "condense_0021_.*.nc"
filename_out1 = "figs/env_a_0021_15.pdf" 

print dir_name1
print filename_in1, filename_out1

make_plot(dir_name1, dir_name2, filename_in1, filename_out1)
#make_plot(dir_name1, filename_in1, filename_out1)
