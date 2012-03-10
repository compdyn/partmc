#!/usr/bin/env python

import sys, re

broken_line_re = re.compile("^(.*[^ ]) *& *")

for filename in sys.argv:
    f = open(filename)
    linenum = 0
    prev_line_broken = False
    prev_line_length = 0
    for line in f:
	linenum += 1
        if prev_line_broken:
            if prev_line_length + len(line.strip()) <= 79:
                print "%s:%d" % (filename, linenum)
        m = broken_line_re.search(line)
        if m:
            prev_line_broken = True
            prev_line_length = len(m.group(1))
        else:
            prev_line_broken = False
    f.close()
