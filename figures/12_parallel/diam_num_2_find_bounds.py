#!/usr/bin/env python
# Copyright (C) 2007-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import numpy as np
import scipy.io
import config
import config_filelist

data_base_dir = "data"
data_type = "diam_num"

if __name__ == "__main__":
    value_min = None
    value_max = None
    for run in config_filelist.runs:
        data_dir = os.path.join(data_base_dir, run["name"])
        for loop in run["loops"]:
            for index in loop["indices"]:
                data_name = "%s_%04d_%08d" % (data_type, loop["num"], index["num"])
                print run["name"] + " " + data_name
                data_filename = os.path.join(data_dir, data_name + ".txt")
                value = np.loadtxt(data_filename)
                mask = np.ma.make_mask(value <= 0.0)
                value = np.ma.array(value, mask=mask)
                print value.min(), value.max()
                print value_min, value_max
                if value_min is None:
                    value_min = value.min()
                else:
                    value_min = min(value_min, value.min())
                if value_max is None:
                    value_max = value.max()
                else:
                    value_max = max(value_max, value.max())
    print value_min, value_max
