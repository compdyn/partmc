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
    '    rgas = const%univ_gas_const * 1d3 / const%air_std_press ! J/mole/K to atmos/(mol/liter)/K',
    '! Copyright (C) 2007-2012, 2016, 2017, 2018, 2021 Nicole Riemer and Matthew West',
    '!!     href="http://dx.doi.org/10.5194/acp-14-5327-2014">10.5194/acp-14-5327-2014</a>.',
    '!!     href="http://dx.doi.org/10.5194/acp-13-11423-2013">10.5194/acp-13-11423-2013</a>.',
    '!!     href="http://dx.doi.org/10.1080/02786826.2020.1804523">10.1080/02786826.2020.1804523</a>.',
    '!!     href="http://dx.doi.org/10.1080/02786826.2019.1661959">10.1080/02786826.2019.1661959</a>',
    '!!     href="http://dx.doi.org/10.5194/gmd-10-4057-2017">10.5194/gmd-10-4057-2017</a>.',
    '!!     href="http://dx.doi.org/10.1080/02786826.2017.1311988">10.1080/02786826.2017.1311988</a>.',
    '!!     href="http://dx.doi.org/10.5194/acp-17-7445-2017">10.5194/acp-17-7445-2017</a>.',
    '!!     href="http://dx.doi.org/10.1016/j.jcp.2016.06.029">10.1016/j.jcp.2016.06.029</a>.',
    '!!     href="http://dx.doi.org/10.5194/acp-14-6289-2014">10.5194/acp-14-6289-2014</a>.',
    ]

for filename in sys.argv[1:]:
    f = open(filename)
    linenum = 0
    for line in f:
        linenum += 1
        if len(line) > 80:
            stripped_line = line.rstrip()
            if not long_comment_re.search(line) and stripped_line not in known_long_lines:
                print("%s:%d" % (filename, linenum))
    f.close()
