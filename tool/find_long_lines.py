#!/usr/bin/env python

import sys, re

long_comment_re = re.compile("^ *!>")

known_long_lines = [
    '!!     href="http://dx.doi.org/10.1016/j.jcp.2011.07.027">10.1016/j.jcp.2011.07.027</a>.',
    '!!     href="http://dx.doi.org/10.1016/j.jaerosci.2009.08.009">10.1016/j.jaerosci.2009.08.009</a>',
    '!!     href="http://dx.doi.org/10.1088/1742-6596/125/1/012020">10.1088/1742-6596/125/1/012020</a>',
    '!> \page input_format_particle Input File Format: Particle-Resolved Simulation',
    '    !! t_progress 600                  # progress printing interval (0 disables) (s)',
    '  !! &lt;name&gt; &lt;whitespace&gt; &lt;data1&gt; &lt;whitespace&gt; &lt;data2&gt; ... # optional comment',
    '! Fortran 95 getopt() and getopt_long(), similar to those in standard C library.',
    '  ! This is needed because Fortran standard allows but doesn\'t *require* short-circuited',
    '                print \'(a,a,a)\', "Error: option \'", trim(arg), "\' requires an argument"',
    ]

for filename in sys.argv[1:]:
    f = open(filename)
    linenum = 0
    for line in f:
	linenum += 1
	if len(line) > 80:
            stripped_line = line.rstrip()
            if not long_comment_re.search(line) and stripped_line not in known_long_lines:
                print "%s:%d" % (filename, linenum)
    f.close()
