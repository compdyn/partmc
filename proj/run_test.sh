#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

#echo "../build/partmc run_mc.spec"
#../build/partmc run_mc.spec

#../build/extract_aero_size_num 1e-9 1e-6 100 out/brownian_mc_0001_ out/brownian_mc_size_num.txt

for f in out/brownian_mc_????_00000001.nc ; do
    f1=${f/_00000001.nc/}
    f2=${f1/_mc_/_mc_size_num_}.txt
    echo "../build/extract_aero_size_num 1e-10 1e-4 220 ${f1}_ $f2"
    ../build/extract_aero_size_num 1e-10 1e-4 220 ${f1}_ $f2
done

for f in out/brownian_mc_????_00000001.nc ; do
    f1=${f/_00000001.nc/}
    f2=${f1/_mc_/_mc_size_mass_}.txt
    echo "../build/extract_aero_size_num 1e-10 1e-4 220 ${f1}_ $f2"
    ../build/extract_aero_size_mass 1e-10 1e-4 220 ${f1}_ $f2
done
exit $?
