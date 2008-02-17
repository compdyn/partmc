#!/bin/sh

cat README
sleep 1

echo ../../src/partmc run_mc.spec
../../src/partmc run_mc.spec

echo ../../src/partmc -p process.dat out/mosaic_post.nc out/mosaic_state_0001_????????.d
../../src/partmc -p process.dat out/mosaic_post_0001 out/mosaic_state_0001_????????.d

echo ./plot.py
./plot.py
echo Now view out/mosaic.pdf
