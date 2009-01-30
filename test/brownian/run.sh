#!/bin/sh

cat README
sleep 1

echo ../../build/partmc run_mc.spec
../../build/partmc run_mc.spec
echo ../../build/partmc run_sect.spec
../../build/partmc run_sect.spec

echo "../../tool/average_netcdf.py -o out/brown_mc_avg.nc out/brown_mc_????.nc"
../../tool/average_netcdf.py -o out/brown_mc_avg.nc out/brown_mc_????.nc

echo ./plot.py
./plot.py
echo "Now view out/brown_*.pdf"
