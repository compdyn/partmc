#!/usr/bin/env python
import subprocess
import os

if not os.path.exists("out"):
    os.mkdir("out")

str_exec = "../../build/partmc"

for hour in range(1, 50):
    print("hour = ", hour)
    
    spec_file_ref = "spec/cond_%02d_ref.spec" % hour
    spec_file_comp = "spec/cond_%02d_comp.spec" % hour 
    spec_file_size = "spec/cond_%02d_size.spec" % hour
    spec_file_both = "spec/cond_%02d_both.spec" % hour

    command_1 = [str_exec, spec_file_ref]
    command_2 = [str_exec, spec_file_comp]
    command_3 = [str_exec, spec_file_size]
    command_4 = [str_exec, spec_file_both]

    print(" ".join(command_1))
    subprocess.check_call(command_1)

    print(" ".join(command_2))
    subprocess.check_call(command_2)

    print(" ".join(command_3))
    subprocess.check_call(command_3)

    print(" ".join(command_4))
    subprocess.check_call(command_4)
