#!/usr/bin/env python

for h in [1, 7, 15, 24]:
    print "%d hour => file %d with start time %d and temp times %d %d" \
          % (h, 6 * h + 1, 3600 * h + 21600, 3600 * h, 3600 * h + 1200)
