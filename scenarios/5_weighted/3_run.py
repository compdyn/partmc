#!/usr/bin/env python
import subprocess

str_exec = "../../build/partmc"

for counter in ["10K", "100K"]:
    print "counter = ", counter
    
    spec_file_flat = "spec/urban_plume_wc_%s_flat.spec" % (counter)
    spec_file_wei1 = "spec/urban_plume_wc_%s_wei-1.spec" % (counter) 
    spec_file_wei2 = "spec/urban_plume_wc_%s_wei-2.spec" % (counter)
    spec_file_wei3 = "spec/urban_plume_wc_%s_wei-3.spec" % (counter)

    command_1 = [str_exec, spec_file_flat]
    command_2 = [str_exec, spec_file_wei1]
    command_3 = [str_exec, spec_file_wei2]
    command_4 = [str_exec, spec_file_wei3]

    print command_1
    subprocess.check_call(command_1)

    print command_2
    subprocess.check_call(command_2)

    print command_3
    subprocess.check_call(command_3)

    print command_4
    subprocess.check_call(command_4)
