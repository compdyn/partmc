#!/usr/bin/env python

n_part = "32"

def process_times(coag_method, output_type):
    times_name = "times_%s_%s_%s.txt" % (n_part, coag_method, output_type)
    times_f = open(times_name, "w")
    #for n in [0, 1, 2, 4, 6, 8, 12, 16]:
    for n in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]:
        try:
            job_filename = "job.urban_plume_%s_%s_%s_%04d.out" \
                       % (n_part, coag_method, output_type, n)
            print job_filename
            job_f = open(job_filename)
            line_p = ""
            line_pp = ""
            for line in job_f:
                if line.startswith("Run finished: timing info above"):
                    cputime = float(line_pp.split()[5])
                    walltime_split = line_p.split()[2].split(":")
                    if len(walltime_split) == 2:
                        walltime = float(walltime_split[0]) * 60 \
                                   + float(walltime_split[1])
                    elif len(walltime_split) == 3:
                        walltime = float(walltime_split[0]) * 3600 \
                                   + float(walltime_split[1]) * 60 \
                                   + float(walltime_split[2])
                    else:
                        raise Exception("unknown walltime format: %s" % str(walltime_split))
                    times_f.write("%d %f %f\n" % (n, cputime, walltime))
                line_pp = line_p
                line_p = line
            job_f.close()
        except Exception, e:
            print "Error: " + job_filename + ": " + str(e)
    times_f.close()

#for coag_method in ["local1", "local2", "local3", "collect",
#                    "central", "dist"]:
for coag_method in ["dist", "local1", "local2", "local3", "none"]:
    if coag_method == "dist":
        #for output_type in ["single", "central", "dist", "none"]:
        for output_type in ["dist"]:
            process_times(coag_method, output_type)
    else:
        output_type = "dist"
        process_times(coag_method, output_type)
