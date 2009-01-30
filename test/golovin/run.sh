#!/bin/bash

cat README
sleep 1

echo ../../build/partmc run_mc.spec
../../build/partmc run_mc.spec
echo ../../build/partmc run_exact.spec
../../build/partmc run_exact.spec

echo "../../tool/average_netcdf.py -o out/golovin_mc_avg.nc out/golovin_mc_????.nc"
../../tool/average_netcdf.py -o out/golovin_mc_avg.nc out/golovin_mc_????.nc
echo ./plot.py
./plot.py
echo Now view out/golovin_num.pdf and out/golovin_vol.pdf
