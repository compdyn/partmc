#!/usr/bin/env python

def process_times(coag_method, output_type):
    times_name = "times_%s_%s.txt" % (coag_method, output_type)
    times_f = open(times_name, "w")
    for n in [0, 1, 2, 4, 6, 8, 12, 16]:
        try:
            job_filename = "job.urban_plume_%s_%s_%02d.out" \
                       % (coag_method, output_type, n)
            print job_filename
            job_f = open(job_filename)
            line_p = ""
            line_pp = ""
            for line in job_f:
                if line.startswith("Run finished: timing info above"):
                    cputime = float(line_pp.split()[5])
                    walltime_min_sec = line_p.split()[2].split(":")
                    walltime = float(walltime_min_sec[0]) * 60 + float(walltime_min_sec[1])
                    times_f.write("%d %f %f\n" % (n, cputime, walltime))
                line_pp = line_p
                line_p = line
            job_f.close()
        except Exception, e:
            print job_filename + ": error"
    times_f.close()

for coag_method in ["local1", "local2", "local3", "collect",
                    "central", "dist"]:
    if coag_method == "dist":
        for output_type in ["single", "central", "dist", "none"]:
            process_times(coag_method, output_type)
    else:
        output_type = "dist"
        process_times(coag_method, output_type)
