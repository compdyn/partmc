#!/usr/bin/env python

import os
import numpy as np
import shutil

if not os.path.exists("spec"):
    os.mkdir("spec")
else:
    shutil.rmtree("spec")
    os.mkdir("spec")

case = 0
for prefactor in np.arange(0.005,0.055,0.005):
    filename_in = "barrel_template.spec"
    for exponent in np.arange(0.1,0.31,0.01): 
        print "prefactor = ", prefactor, "exponent = ", exponent
        case += 1
        filename_out = "spec/barrel_wc_case_%04d.spec" % (case)
        f_in = open(filename_in, 'r')
        f_out = open(filename_out, 'w')

        for line in f_in:
            line = line.replace('%%prefactor%%', str(prefactor))
            line = line.replace('%%exponent%%', str(exponent))
            f_out.write(line)

        f_in.close()
        f_out.close()
