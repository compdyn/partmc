#!/usr/bin/env python
import subprocess
import os
import math
from numpy import *

if not os.path.exists("spec_with_bounds"):
    os.mkdir("spec_with_bounds")

bound_list = ['0925_low_size_high_conc','0925_low_size_low_conc',\
              '0925_high_size_high_conc','0925_high_size_low_conc','0925']

if not os.path.exists("out_with_bounds"):
       os.mkdir("out_with_bounds")

prefactor = float(raw_input("Enter prefactor:"))
exponent = float(raw_input("Enter exponent:"))
filename_in = "barrel_template.spec"
for prefix in bound_list:
    filename_out = "spec_with_bounds/barrel_wc_"+prefix+".spec"
    f_in = open(filename_in, 'r')
    f_out = open(filename_out, 'w')
    for line in f_in:
        line = line.replace('%%OUTPUT_PREFIX%%', 'out_with_bounds/'+prefix)
        line = line.replace('%%prefactor%%', str(prefactor))
        line = line.replace('%%exponent%%', str(exponent))
        f_out.write(line)

    f_in.close()
    f_out.close()

str_exec = "../../build/partmc"
str_extr_size = "../../build/extract_aero_size"
str_extr_time = "../../build/extract_aero_time"

for prefix in bound_list:
    command_0 = ["cp", "aero_init_size_dist_"+prefix+".dat", "aero_init_size_dist.dat"]
    subprocess.check_call(command_0)
    ncfile_prefix = "out_with_bounds/"+prefix+"_0001"
    spec_file = "spec_with_bounds/barrel_wc_"+prefix+".spec"
    command_1 = [str_exec,spec_file]
    command_2 = [str_extr_size,"--num","--dmin","1e-8","--dmax","1e-6","--nbin","100",ncfile_prefix]
    command_3 = [str_extr_size,"--mass","--dmin","1e-8","--dmax","1e-6","--nbin","100",ncfile_prefix]
    command_4 = [str_extr_time,ncfile_prefix]

    print command_1
    subprocess.check_call(command_1)
    print command_2
    subprocess.check_call(command_2)
    print command_3
    subprocess.check_call(command_3)
    print command_4
    subprocess.check_call(command_4)
