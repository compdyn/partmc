#!/usr/bin/env python
import subprocess

str_exec_size = "../../build/bin_average_size"
str_exec_comp = "../../build/bin_average_comp"

bin_start = "1e-10"
bin_end = "1e-5"
bin_n ="25" 
flag_dry_wet = "dry"
flag_position = "average"
flag_measurement = "number"

for hour in range(1, 50):
    print("hour: %d" % hour)
    
    file_in = "start/urban_plume2_0001_000000%02d.nc" % hour
    file_in_comp = "start/urban_plume2_comp_0001_000000%02d.nc" % hour
    file_out_size = "start/urban_plume2_size"
    file_out_comp = "start/urban_plume2_comp"
    file_out_both = "start/urban_plume2_both"

    command_1 = [str_exec_size, bin_start, bin_end, bin_n, flag_dry_wet, flag_position, flag_measurement, file_in, file_out_size]
    command_2 = [str_exec_comp, bin_start, bin_end, bin_n, flag_dry_wet, file_in, file_out_comp]
    command_3 = [str_exec_size, bin_start, bin_end, bin_n, flag_dry_wet, flag_position, flag_measurement, file_in_comp, file_out_both]
    print(" ".join(command_1))
    subprocess.check_call(command_1)
    
    print(" ".join(command_2))
    subprocess.check_call(command_2)

    print(" ".join(command_3))
    subprocess.check_call(command_3)
