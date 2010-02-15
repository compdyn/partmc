#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

maximum_ss = np.zeros([4,49])
 
def make_plot(in_dir, in_file_pattern, out_filename, title, max_ss_i, max_ss_j):
    print in_dir, in_file_pattern
    
    env_state_history = partmc.read_history(partmc.env_state_t, in_dir, in_file_pattern)
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

for counter in range(1,41):
    print "counter = ", counter
    
    dir_name = "../../scenarios/3_condense/out"
    filename_in1 = "cond_%02d_ref_0001_.*.nc" % counter
    filename_in2 = "cond_%02d_comp_0001_.*.nc" % counter
    filename_in3 = "cond_%02d_size_0001_.*.nc" % counter
    filename_in4 = "cond_%02d_both_0001_.*.nc" % counter

    filename_out1 = "figs/env_wc_ref_%02d.pdf" % (counter-1)
    filename_out2 = "figs/env_wc_comp_%02d.pdf" % (counter-1)
    filename_out3 = "figs/env_wc_size_%02d.pdf" % (counter-1)
    filename_out4 = "figs/env_wc_both_%02d.pdf" % (counter-1)

    title = " %02d hours" % (counter-1)

    print dir_name, title
    print filename_in1, filename_out1
    print filename_in2, filename_out2
    print filename_in3, filename_out3
    print filename_in4, filename_out4

    make_plot(dir_name,filename_in1, filename_out1, title, 0, counter-1)
    make_plot(dir_name,filename_in2, filename_out2, title, 1, counter-1)
    make_plot(dir_name,filename_in3, filename_out3, title, 2, counter-1)
    make_plot(dir_name,filename_in4, filename_out4, title, 3, counter-1)

np.savetxt("data/maximum_ss.txt", maximum_ss)
