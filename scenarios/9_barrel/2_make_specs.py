#!/usr/bin/env python

import os
import numpy as np

if not os.path.exists("spec"):
    os.mkdir("spec")

run = 0
for prefactor in np.arange(0.005,0.051,0.001):
    filename_in = "barrel_template.spec"
    for exponent in np.arange(0.2,0.31,0.01): 
        print "prefactor = ", prefactor, "exponent = ", exponent
        run = run + 1
        filename_out = "spec/barrel_wc_run_%04d.spec" % (run)
        f_in = open(filename_in, 'r')
        f_out = open(filename_out, 'w')

        for line in f_in:
            line = line.replace('%%prefactor%%', str(prefactor))
            line = line.replace('%%exponent%%', str(exponent))
            f_out.write(line)

        f_in.close()
        f_out.close()



