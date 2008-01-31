#!/usr/bin/env python

import sys

for filename in sys.argv:
    f = open(filename)
    linenum = 0
    for line in f:
	linenum += 1
	if len(line) > 80:
	    print "%s:%d" % (filename, linenum)
    f.close()
