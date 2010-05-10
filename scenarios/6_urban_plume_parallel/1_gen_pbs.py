#!/usr/bin/env python

import math

def sub_file(in_filename, out_filename, subs):
    f_in = open(filename_in, 'r')
    f_out = open(filename_out, 'w')
    for line in f_in:
        for (template, replacement) in subs.iteritems():
            line = line.replace(template, replacement)
        f_out.write(line)
    f_in.close()
    f_out.close()

######################################################################
    
spec_template_filename = "urban_plume_template.spec"

for coag_method in ["local1", "local2", "local3", "collect",
                    "central", "dist"]:
    for output_method in ["single", "central", "dist", "none"]:
        if coag_method.startswith("local"):
            if coag_method == "local1":
                mix_prob = 0.01
            if coag_method == "local2":
                mix_prob = 0.1
            if coag_method == "local3":
                mix_prob = 0.5
            del_t = 60
            mix_timescale = "%.2f" \
                            % (- del_t / math.log(1 - mix_prob))
        else:
            mix_timescale = "0"
        out_filename = "urban_plume_%s_%s.spec" \
                       % (coag_method, output_method)
        sub_file(spec_template_filename, out_filename,
                 {"%%{{OUTPUT_TYPE}}%%": output_type,
                  "%%{{MIX_TIMESCALE}}%%": mix_timescale,
                  "%%{{COAG_METHOD}}%%": coag_method,
                  })

######################################################################

pbs_template_filename = "run_mvapich-intel_template.pbs"
max_walltime = "01:00:00" # hh:mm:ss

for coag_method in ["local1", "local2", "local3", "collect",
                    "central", "dist"]:
    for output_method in ["single", "central", "dist", "none"]:
        spec_file_name = "urban_plume_with_coag_%s_%s.spec" \
                         % (coag_method, output_method)
        for n in [0, 1, 2, 4, 6, 8, 12, 16]:
            job_name = "urban_plume_%s_%s_%02d" \
                       % (coag_method, output_method, n)
            pbs_file_name = "run_%s.pbs" % job_name
            if n == 0:
                num_nodes = "1"
                procs_per_node = "1"
            else:
                num_nodes = str(n)
                procs_per_node = "8"
            sub_file(pbs_template_filename, pbs_file_name,
                     {"%%{{MAX_WALLTIME}}%%": max_walltime,
                      "%%{{NUM_NODES}}%%": num_nodes,
                      "%%{{PROCS_PER_NODE}}%%": procs_per_node,
                      "%%{{SPEC_FILE_NAME}}%%": spec_file_name,
                      "%%{{JOB_NAME}}%%": job_name,
                      })
