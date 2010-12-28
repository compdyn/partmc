#!/usr/bin/env python

import os

if not os.path.exists("spec"):
    os.mkdir("spec")
if not os.path.exists("temp"):
    os.mkdir("temp")

for counter in range(1, 50):
    print "counter = ", counter
    filename_in = "cond_template.spec"
    for run in ["ref", "comp", "size", "both"]: 
        filename_out = "spec/cond_%02d_%s.spec" % (counter, run)
        print "filename_out ", filename_out
        f_in = open(filename_in, 'r')
        f_out = open(filename_out, 'w')

        for line in f_in:
            line = line.replace('%%OUTPUT_PREFIX%%', 'out/cond_%02d_%s' % (counter, run))
            if run == "ref":
                line = line.replace('%%RESTART_FILE%%',  'start/urban_plume_wc_0001_000000%02d.nc' % counter )
            else:
                line = line.replace('%%RESTART_FILE%%',  'start/urban_plume_%s_wc_0001_000000%02d.nc' % (run, counter) )
            line = line.replace('%%TEMP_PROFILE%%',  'temp/temp_%02d.dat' % counter)
            f_out.write(line)

        f_in.close()
        f_out.close()

    filename_out_temp = "temp/temp_%02d.dat" % counter
    print "filename_out_temp", filename_out_temp
    f_out = open(filename_out_temp, 'w')

    f_out.write("# time (s)\n")
    f_out.write("# temp (K)\n")
    f_out.write("time  %.2f %.2f\n" % ((counter-1) * 3600.0, (counter-1) * 3600.0 + 1200))
    f_out.write("temp  290   280\n")   
    f_out.close()

