#!/usr/bin/env python

import os

if not os.path.exists("spec"):
    os.mkdir("spec")

for counter in ["10K", "100K"]:
    print "counter = ", counter
    filename_in = "weighted_template.spec"
    for run in ["flat", "wei-1", "wei-2", "wei-3"]: 
        filename_out = "spec/urban_plume_wc_%s_%s.spec" % (counter, run)
        print "filename_out ", filename_out
        f_in = open(filename_in, 'r')
        f_out = open(filename_out, 'w')

        for line in f_in:
            line = line.replace('%%OUTPUT_PREFIX%%', 'out/urban_plume_wc_%s_%s' % (counter, run))
            if run == "flat":
                line = line.replace('%%WEIGHTING_FUNC%%',  'none               ')
                if line.find('%%EXPONENT%%') > -1:
                    continue
                if line.find('%%REF_RADIUS%%') > -1:
                    continue
            else:
                line = line.replace('%%WEIGHTING_FUNC%%',  'power              ')
                line = line.replace('%%REF_RADIUS%%',  '1e-7              ')
                if run == "wei-1":
                    line = line.replace('%%EXPONENT%%',  '-1              ')
                elif run == "wei-2":
                    line = line.replace('%%EXPONENT%%',  '-2              ') 
                elif run == "wei-3":
                    line = line.replace('%%EXPONENT%%',  '-3              ')
            f_out.write(line)

        f_in.close()
        f_out.close()



