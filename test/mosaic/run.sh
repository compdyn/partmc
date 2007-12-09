#!/bin/sh

cat <<ENDINFO

MOSAIC Test-case
----------------

This tests the interface to the MOSAIC chemistry code.

ENDINFO
sleep 1

echo ../../src/partmc run_mc.spec
../../src/partmc run_mc.spec

echo ./plot.py
./plot.py
echo Now view out/mosaic.pdf
