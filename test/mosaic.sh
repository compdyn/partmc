#!/bin/sh

cat <<ENDINFO

MOSAIC Test-case
----------------

This tests the interface to the MOSAIC chemistry code.

ENDINFO
sleep 1

echo ../src/partmc mosaic.spec
../src/partmc mosaic.spec

echo ./mosaic_plot.py
./mosaic_plot.py
echo Now view out/mosaic.pdf
