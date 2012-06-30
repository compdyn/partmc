#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_dir, in_file_pattern, out_filename):
    print in_dir, in_file_pattern
    
    env_state_history = partmc.read_history(partmc.env_state_t, in_dir, in_file_pattern)
    time = [env_state_history[i][0] for i in range(len(env_state_history))]
    rh = [env_state_history[i][1].relative_humidity for i in range(len(env_state_history))]
    temp = [env_state_history[i][1].temperature for i in range(len(env_state_history))]

    print (max(rh) - 1)*100.
    
    plt.clf()
    plt.plot(time,rh,'r')
    ax1 = plt.gca()
    ax1.set_ylim(0.9,1.003)
    ax2 = plt.twinx()
    ax2.plot(time,temp)
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/8_condense_p/out/"

filename_in1 = "condense_0001_.*.nc"
filename_out1 = "figs/env_ne.pdf" 

print dir_name
print filename_in1, filename_out1

make_plot(dir_name,filename_in1, filename_out1)

