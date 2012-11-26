#!/usr/bin/env python
import subprocess
import os
import math
from numpy import *

def read_values_from_spec_file(filename_in, str_wanted):
    A = []
    f_in = open(filename_in, 'r')
    for line in f_in:
        s =  line.rsplit(' ')
        if (s[0] == str_wanted):
           A.append(s[1]) # add the value to the list
    f_in.close()
    # Convert strings to an array of values (in this case floats)
    val = zeros(len(A))
    for ii in range(0,len(A)):
        val[ii] =  float(A[ii])
    return val

if not os.path.exists("out"):
    os.mkdir("out")

str_exec = "../../build/partmc"
str_extr_num = "../../build/extract_aero_size"
str_extr_mass = "../../build/extract_aero_size"

# Open file to store the errors
filename_out_num = "rel_err_num.dat"
f_out_num = open(filename_out_num, 'w')
f_out_num.write("# Colume 1: prefactor\n")
f_out_num.write("# Colume 2: exponent\n")
f_out_num.write("# Colume j+2: relative error of number distribution at time(j)\n")
filename_out_mass = "rel_err_mass.dat"
f_out_mass = open(filename_out_mass, 'w')
f_out_mass.write("# Colume 1: prefactor\n")
f_out_mass.write("# Colume 2: exponent\n")
f_out_mass.write("# Colume j+2: relative error of mass distribution at time(j)\n")
filename_out_total = "rel_err_total.dat"
f_out_total = open(filename_out_total, 'w')
f_out_total.write("# Colume 1: prefactor\n")
f_out_total.write("# Colume 2: exponent\n")
f_out_total.write("# Colume j+2: total relative error at time(j)\n")
path="spec/"
dirList=os.listdir(path)
run = 0
for fname in dirList:

    run += 1
    spec_file = "spec/barrel_wc_run_%04d.spec" % (run)

    # Read prefactor and expoenent values from spec file
    prefactor = read_values_from_spec_file(spec_file, 'prefactor_BL')
    exponent = read_values_from_spec_file(spec_file, 'exponent_BL')

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

    f_out_num.write("%.3f      %.3f     " % (prefactor, exponent))
    f_out_mass.write("%.3f      %.3f     " % (prefactor, exponent))
    f_out_total.write("%.3f      %.3f     " % (prefactor, exponent))
    for col in range(1,data1.shape[1]):
        data1_1d = data1[:,col]
        data2_1d = data2[:,col]
        data3_1d = data3[:,col]
        data4_1d = data4[:,col]
        # extract_aero_size gives d*_dlnDp, neet to convert to d*_dlogDp
        data1_1d *= math.log(10)
        data3_1d *= math.log(10)
        # calculate and store relative error
        diff = data2_1d - data1_1d
        rel_err_num = sqrt(sum(diff**2)) / sqrt(sum(data2_1d**2))
        diff = data4_1d - data3_1d
        rel_err_mass = sqrt(sum(diff**2)) / sqrt(sum(data4_1d**2))
        rel_err_total = sqrt(rel_err_num**2 + rel_err_mass**2)
        if col != data1.shape[1]-1:
           f_out_num.write("%.3f     " % (rel_err_num))
           f_out_mass.write("%.3f     " % (rel_err_mass))
           f_out_total.write("%.3f     " % (rel_err_total))
        else:
           f_out_num.write("%.3f\n" % (rel_err_num))
           f_out_mass.write("%.3f\n" % (rel_err_mass))
           f_out_total.write("%.3f\n" % (rel_err_total))

f_out_num.close()
f_out_mass.close()
f_out_total.close()
