#!/usr/bin/env python

import subprocess
import os
import config

list_exec = ["mpirun", "-v", "-np", "10", "../../../../build/partmc"]

for run in config.all_runs():
    print run["name"]
    dirname = os.path.join(config.run_dirname, run["name"])
    subprocess.check_call(list_exec + ["run_sect.spec"], cwd=dirname)
    subprocess.check_call(list_exec + ["run_part.spec"], cwd=dirname)
