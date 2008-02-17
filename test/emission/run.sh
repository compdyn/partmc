#!/bin/sh

cat README
sleep 1

echo ../../src/partmc run_mc.spec
../../src/partmc run_mc.spec
echo ../../src/partmc run_sect.spec
../../src/partmc run_sect.spec
echo ../../src/partmc run_exact.spec
../../src/partmc run_exact.spec

echo ./plot.py
./plot.py
echo Now view out/emission_dist.pdf and out/emission_history.pdf
