#!/usr/bin/env python

import sys, re

for filename in sys.argv:
    f = open(filename)
    linenum = 0
    for line in f:
	linenum += 1
	if len(line.rstrip('\r\n')) > len(line.rstrip()):
            print "%s:%d" % (filename, linenum)
    f.close()
