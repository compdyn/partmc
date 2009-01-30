#!/bin/sh

cat README
sleep 1

echo ../../build/partmc run_mc.spec
../../build/partmc run_mc.spec
echo ../../build/partmc run_sect.spec
../../build/partmc run_sect.spec

echo ./plot.py
./plot.py
echo "Now view out/sedi_*.pdf"
