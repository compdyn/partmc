#!/usr/bin/env python
# Copyright (C) 2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import sys
sys.path.append("../tool")
from pmc_data_nc import *

n_bin = 450
ss_active_axis = pmc_linear_axis(0.001, 0.01, n_bin)
