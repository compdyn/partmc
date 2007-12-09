#!/bin/sh

cat <<ENDINFO

MOSAIC Test-case
----------------

This tests the interface to the MOSAIC chemistry code. The number of
aerosol particles is fixed at 3 and there is no coagulation. Only
chemistry and aerosol-gas transfers occur.

ENDINFO
sleep 1

echo ../../src/partmc run_mc.spec
../../src/partmc run_mc.spec

echo ./plot.py
./plot.py
echo Now view out/mosaic.pdf
