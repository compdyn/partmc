#!/usr/bin/env python
import subprocess
import os
import math
from numpy import *

if not os.path.exists("out"):
    os.mkdir("out")

str_exec = "../../build/partmc"
str_extr_num = "../../build/extract_aero_size"
str_extr_mass = "../../build/extract_aero_size"

# Open file to store the errors
filename_out = "rel_err.dat"
f_out = open(filename_out, 'w')
f_out.write("run_id  rel_err\n")
for run in range(1,507):
    print "run = ", run

    spec_file = "spec/barrel_wc_run_%04d.spec" % (run)

    command_1 = [str_exec,spec_file]
    command_2 = [str_extr_num,"--num","--dmin","1e-8","--dmax","1e-6","--nbin","100","out/barrel_wc_0001"]
    command_3 = [str_extr_mass,"--mass","--dmin","1e-8","--dmax","1e-6","--nbin","100","out/barrel_wc_0001"]

    print command_1
    subprocess.check_call(command_1)

    print command_2
    subprocess.check_call(command_2)
    data1 = loadtxt("out/barrel_wc_0001_aero_size_num.txt")
    data2 = loadtxt("ref_aero_size_num_regrid.txt")

    print command_3
    subprocess.check_call(command_3)
    data3 = loadtxt("out/barrel_wc_0001_aero_size_mass.txt")
    data4 = loadtxt("ref_aero_size_mass_regrid.txt")

    max_rel_err = 0.
    for col in range(0,data1.shape[1]):
        data1_1d = data1[:,col]
        data2_1d = data2[:,col]
        data3_1d = data3[:,col]
        data4_1d = data4[:,col]
        # extract_aero_size gives d*_dlnDp, neet to convert to d*_dlogDp
        if (col != 0):
            data1_1d *= math.log(10)
            data3_1d *= math.log(10)
        # calculate relative error, store the max value
        diff = data2_1d - data1_1d
        rel_err_num = sqrt(sum(diff**2)) / sqrt(sum(data2_1d**2))
        diff = data4_1d - data3_1d
        rel_err_mass = sqrt(sum(diff**2)) / sqrt(sum(data4_1d**2))
        max_rel_err = max(max_rel_err, max(rel_err_num, rel_err_mass))

    f_out.write("%04d    %.3f\n" % (run, max_rel_err))
    if (run == 1):
        min_err = max_rel_err
    else:    
        min_err = min(min_err, max_rel_err)
    f_out.write("Minimum error is: %.3f, occurred at run %04d\n" % (min_err, run))

f_out.write("End of Run: minimum error is: %.3f, occurred at run %04d\n" % (min_err, run))
f_out.close()
