#!/usr/bin/env python

import os
import shutil
import config

def copy_template_file(filename_in, filename_out, run):
    f_in = open(filename_in, 'r')
    f_out = open(filename_out, 'w')
    for line in f_in:
        line = line.replace('%%N_PART%%', str(int(run["n_part"]) * config.n_nodes * config.n_procs_per_node))
        if run["weight_type"] == "power":
            line = line.replace('%%WEIGHT_TYPE%%', run["weight_type"] + "\nexponent " + run["exponent"])
        else:
            line = line.replace('%%WEIGHT_TYPE%%', run["weight_type"])
        f_out.write(line)
    f_out.close()
    f_in.close()

if not os.path.exists(config.run_dirname):
    os.mkdir(config.run_dirname)
for run in config.all_runs():
    print run["name"]
    dirname = os.path.join(config.run_dirname, run["name"])
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    for f in [
        "aero_back.dat",
        "aero_back_dist.dat",
        "aero_data.dat",
        "aero_emit.dat",
        "aero_emit_dist.dat",
        "aero_init_comp.dat",
        "aero_init_dist.dat",
        "gas_back.dat",
        "gas_data.dat",
        "gas_emit.dat",
        "gas_init.dat",
        "height.dat",
        "run_part.spec",
        "run_sect.spec",
        "temp.dat",
        ]:
        copy_template_file(os.path.join(config.ref_dirname, f),
                           os.path.join(dirname, f), run)
