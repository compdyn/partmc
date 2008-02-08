#!/bin/sh

cat <<ENDINFO

MOSAIC Test-case
----------------

This tests the interface to the MOSAIC chemistry code. The number of
aerosol particles is fixed at 3 and there is no coagulation. Only
chemistry and aerosol-gas transfers occur. This also tests the full
state output and post-processing capabilities of PartMC. The graph
symbols derived from post-processing should overlay the lines.

ENDINFO
sleep 1

echo ../../src/partmc run_mc.spec
../../src/partmc run_mc.spec

echo ../../src/partmc -p process.dat out/mosaic_post.nc out/mosaic_state_0001_????????.d
../../src/partmc -p process.dat out/mosaic_post.nc out/mosaic_state_0001_????????.d

echo ./plot.py
./plot.py
echo Now view out/mosaic.pdf
