#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
import subprocess
import os

matplotlib.use("PDF")
matplotlib.use('Agg')

import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_dir, in_file_pattern, out_filename):
    print in_dir, in_file_pattern
    
    gas_state_history = partmc.read_history(partmc.gas_state_t, in_dir, in_file_pattern)
    time = [t for [t, gs] in gas_state_history]
    o3 = [gs.mixing_ratio('O3') for [t, gs] in gas_state_history]

    plt.clf()
    plt.plot(time,o3,'r')
    fig.savefig(out_filename)

dir_name = "../../scenarios/5_weighted/out"
filename_in1 = "urban_plume_wc_10K_flat_0001_.*.nc" % counter

filename_out1 = "figs/env_ref_%02d.png" % (counter-1)

make_plot(dir_name,filename_in1, filename_out1)

