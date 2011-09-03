#!/usr/bin/env python

import sys, os, re

error_code_re = re.compile("\(([0-9]{9}),")

error_codes = {}

if len(sys.argv) < 2:
    print("usage: {} <files>".format(os.path.basename(sys.argv[0])))

for filename in sys.argv[1:]:
    with open(filename) as f:
        for (line_number, line) in enumerate(f, start=1):
            m = error_code_re.search(line)
            if m:
                code = m.group(1)
                if code not in error_codes:
                    error_codes[code] = []
                error_codes[code].append("{}:{}".format(filename, line_number))
    
for code in error_codes.keys():
    if len(error_codes[code]) > 1:
        print("{}: {}".format(code, ", ".join(error_codes[code])))
