#!/usr/bin/env python

import math, os

if not os.path.exists("specs"):
    os.mkdir("specs")
if not os.path.exists("pbs"):
    os.mkdir("pbs")

def sub_file(in_filename, out_filename, subs):
    f_in = open(in_filename, 'r')
    f_out = open(out_filename, 'w')
    for line in f_in:
        for (template, replacement) in subs.iteritems():
            line = line.replace(template, replacement)
        f_out.write(line)
    f_in.close()
    f_out.close()

######################################################################
    
spec_template_filename = "urban_plume_template.spec"
n_part = "32"

for coag_method_var in ["none", "local1", "local2", "local3", "collect",
                    "central", "dist"]:
    for output_type_var in ["single", "central", "dist", "none"]:
        if coag_method_var == "none":
            coag_method = "local"
            mix_timescale = "0"
            gas_average = "no"
            env_average = "no"
        elif coag_method_var.startswith("local"):
            coag_method = "local"
            if coag_method_var == "local1":
                mix_prob = 0.01
            if coag_method_var == "local2":
                mix_prob = 0.1
            if coag_method_var == "local3":
                mix_prob = 0.9999
            del_t = 60
            mix_timescale = "%.2f" \
                            % (- del_t / math.log(1 - mix_prob))
            gas_average = "yes"
            env_average = "yes"
        else:
            coag_method = coag_method_var
            mix_timescale = "0"
            gas_average = "yes"
            env_average = "yes"
        if output_type_var == "none":
            t_output = "0"
            output_type = "dist"
        else:
            t_output = "3600"
            output_type = output_type_var
        out_filename = "specs/urban_plume_%s_%s_%s.spec" \
                       % (n_part, coag_method_var, output_type_var)
        sub_file(spec_template_filename, out_filename,
                 {"%%{{OUTPUT_TYPE}}%%": output_type,
                  "%%{{MIX_TIMESCALE}}%%": mix_timescale,
                  "%%{{COAG_METHOD}}%%": coag_method,
                  "%%{{T_OUTPUT}}%%": t_output,
                  "%%{{GAS_AVERAGE}}%%": gas_average,
                  "%%{{ENV_AVERAGE}}%%": env_average,
                  "%%{{N_PART}}%%": n_part,
                  })

######################################################################

pbs_template_filename = "run_mvapich-intel_template.pbs"
max_walltime = "04:00:00" # hh:mm:ss

def make_pbs(coag_method, output_type):
    spec_file_name = "specs/urban_plume_%s_%s_%s.spec" \
                     % (n_part, coag_method, output_type)
    for n in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]:
        job_name = "urban_plume_%s_%s_%s_%04d" \
                   % (n_part, coag_method, output_type, n)
        pbs_file_name = "pbs/run_%s.pbs" % job_name
        if n < 8:
            num_nodes = "1"
            procs_per_node = str(n)
        else:
            num_nodes = str(n / 8)
            procs_per_node = "8"
        sub_file(pbs_template_filename, pbs_file_name,
                 {"%%{{MAX_WALLTIME}}%%": max_walltime,
                  "%%{{NUM_NODES}}%%": num_nodes,
                  "%%{{PROCS_PER_NODE}}%%": procs_per_node,
                  "%%{{SPEC_FILE_NAME}}%%": spec_file_name,
                  "%%{{JOB_NAME}}%%": job_name,
                  })

#for coag_method in ["none", "local1", "local2", "local3", "collect",
#                    "central", "dist"]:
for coag_method in ["none", "local1", "local2", "local3", "dist"]:
    if coag_method == "dist":
        #for output_type in ["single", "central", "dist", "none"]:
        for output_type in ["dist"]:
            make_pbs(coag_method, output_type)
    else:
        output_type = "dist"
        make_pbs(coag_method, output_type)
