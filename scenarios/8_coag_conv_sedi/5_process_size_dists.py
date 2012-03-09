#!/usr/bin/env python

import subprocess
import os
import config

part_num_exec = ["../../../../build/extract_aero_size", "--num", "--dmin", "1e-7", "--dmax", "1", "--nbin", "100"]
part_mass_exec = ["../../../../build/extract_aero_size", "--mass", "--dmin", "1e-7", "--dmax", "1", "--nbin", "100"]
sect_num_exec = ["../../../../build/extract_sectional_aero_size", "--num"]
sect_mass_exec = ["../../../../build/extract_sectional_aero_size", "--mass"]

for run in config.all_runs():
    print run["name"]
    dirname = os.path.join(config.run_dirname, run["name"])
    subprocess.check_call(part_num_exec + ["part_0001_0001"], cwd=dirname)
    subprocess.check_call(part_mass_exec + ["part_0001_0001"], cwd=dirname)
    subprocess.check_call(sect_num_exec + ["sect"], cwd=dirname)
    subprocess.check_call(sect_mass_exec + ["sect"], cwd=dirname)
