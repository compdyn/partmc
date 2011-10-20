#!/usr/bin/env python

import subprocess
import os
import config

str_exec = "../../../../build/partmc"

for run in config.all_runs():
    print run["name"]
    dirname = os.path.join(config.run_dirname, run["name"])
    subprocess.check_call([str_exec, "run_sect.spec"], cwd=dirname)
    subprocess.check_call([str_exec, "run_part.spec"], cwd=dirname)
