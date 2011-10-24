#!/usr/bin/env python

import os
import shutil
import config

def copy_template_file(filename_in, filename_out, job_name):
    f_in = open(filename_in, 'r')
    f_out = open(filename_out, 'w')
    for line in f_in:
        line = line.replace('%%JOB_NAME%%', job_name)
        line = line.replace('%%N_NODES%%', str(config.n_nodes))
        line = line.replace('%%N_PROCS_PER_NODE%%', str(config.n_procs_per_node))
        if job_name == "1k":
            line = line.replace('%%WALL_TIME%%', "03:00:00")
        elif job_name == "10k":
            line = line.replace('%%WALL_TIME%%', "30:00:00")
        elif job_name == "100k":
            line = line.replace('%%WALL_TIME%%', "300:00:00")
        f_out.write(line)
    f_out.close()
    f_in.close()

qsubs = set()
for run in config.all_runs():
    dirname = os.path.join(config.run_dirname, run["name"])
    qsubs.add(run["n_part_name"])

for n_part_name in qsubs:
    print n_part_name
    copy_template_file("template.qsub",
                       os.path.join(config.run_dirname,
                                    "run_" + n_part_name + ".qsub"),
                       n_part_name)
