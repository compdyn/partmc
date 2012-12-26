#!/usr/bin/env python
import subprocess

str_exec_size = "../../build/bin_average_size"
str_exec_comp = "../../build/bin_average_comp"

bin_start = "1e-10"
bin_end = "1e-5"
bin_n ="1" 
flag_dry_wet = "dry"
flag_position = "average"
flag_quantity = "number"

for hour in range(1, 26):
    print "hour = ", hour
    for loop in range (1,11):

        file_in = "out/urban_plume_00%02d_000000%02d.nc" % (loop, hour)
        file_in_comp = "out/urban_plume_comp1_00%02d_000000%02d.nc" % (loop, hour)
        file_out_size = "out/urban_plume_size1"
        file_out_comp = "out/urban_plume_comp1"
        file_out_both = "out/urban_plume_both1"

        command_1 = [str_exec_size, bin_start, bin_end, bin_n, flag_dry_wet, flag_position, flag_quantity, file_in, file_out_size]
        command_2 = [str_exec_comp, bin_start, bin_end, bin_n, flag_dry_wet, file_in, file_out_comp]
        command_3 = [str_exec_size, bin_start, bin_end, bin_n, flag_dry_wet, flag_position, flag_quantity, file_in_comp, file_out_both]
        print command_1
        subprocess.check_call(command_1)

        print command_2
        subprocess.check_call(command_2)

        print command_3
        subprocess.check_call(command_3)
