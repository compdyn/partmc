#!/bin/sh

cat README
sleep 1

echo ../../src/partmc run_mc.spec
../../src/partmc run_mc.spec
echo ../../src/partmc run_sect.spec
../../src/partmc run_sect.spec

echo ./plot.py
./plot.py
echo "Now view out/sedi_*.pdf"
