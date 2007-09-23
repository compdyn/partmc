#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *

color_list = [color.hsb(0/3.0, 1, 1),
	      color.hsb(1/3.0, 1, 1),
	      color.hsb(2/3.0, 1, 1),
	      color.hsb(1/6.0, 1, 1),
	      color.hsb(3/6.0, 1, 1),
	      color.hsb(5/6.0, 1, 1),
	      color.hsb(1/12.0, 1, 1),
	      color.hsb(3/12.0, 1, 1),
	      color.hsb(5/12.0, 1, 1),
	      color.hsb(7/12.0, 1, 1),
	      color.hsb(9/12.0, 1, 1),
	      color.hsb(11/12.0, 1, 1),
	      ]

grid_painter = graph.axis.painter.regular(gridattrs = [style.linestyle.dotted])
