#!/bin/sh

cat README
sleep 1

echo ../../build/partmc run_mc.spec
../../build/partmc run_mc.spec

echo ../../build/partmc -p process.dat out/mosaic_post.nc out/mosaic_state_0001_????????.d
../../build/partmc -p process.dat out/mosaic_post_0001 out/mosaic_state_0001_????????.d

echo ./plot.py
./plot.py
echo Now view out/mosaic.pdf
