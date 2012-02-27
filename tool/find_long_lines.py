#!/usr/bin/env python

import sys, re

long_comment_re = re.compile("^ *!>")

for filename in sys.argv:
    f = open(filename)
    linenum = 0
    for line in f:
	linenum += 1
	if len(line) > 80:
            if not long_comment_re.search(line):
                print "%s:%d" % (filename, linenum)
    f.close()
