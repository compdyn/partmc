#!/bin/bash

set -e
set -v
cd ${0%/*}

add_pound() {
    spec_file=$1
    head=$2
    N=`grep -n "^$head" $spec_file | cut -d: -f1`
    sed -i "$N s/^/#/" $spec_file
}    

sub_pound() {
    spec_file=$1
    head=$2
    N=`grep -n "^#$head" $spec_file | cut -d: -f1`
    sed -i "$N s/^#//" $spec_file
}    

comment_for_restart() {
    spec_file=$1
    add_pound $1 gas_data
    add_pound $1 gas_init
    add_pound $1 aerosol_data
    add_pound $1 do_fractal
    add_pound $1 aerosol_init
    add_pound $1 do_select_weighting
    add_pound $1 weight_type
    
}

uncomment_for_restart() {
    spec_file=$1
    sub_pound $1 gas_data
    sub_pound $1 gas_init
    sub_pound $1 aerosol_data
    sub_pound $1 do_fractal
    sub_pound $1 aerosol_init
    sub_pound $1 do_select_weighting
    sub_pound $1 weight_type
    
}

caseName=$1
chi=$2
t_max=$3
frzDir=/data/nriemer/a/wenhant2/modeling/partmc/scenarios/7_freezing
cd $frzDir

mkdir -p $frzDir/output/$caseName
cp -p $frzDir/run_part.spec $frzDir/output/$caseName
cp -p $frzDir/*.dat $frzDir/output/$caseName

sed -i "/output_prefix /coutput_prefix output/${caseName}/freezing_part # prefix of output files" $frzDir/run_part.spec

sed -i "/t_max /ct_max 0                            # total simulation time (s)" $frzDir/run_part.spec
../../build/partmc run_part.spec
$frzDir/aero_init_tools/run_setchi.py $chi $frzDir/output/${caseName}/freezing_part_0001_00000001.nc $frzDir/aero_init_comp.dat

sed -i "/restart /crestart yes                     # whether to restart from saved state (yes/no)" $frzDir/run_part.spec
sub_pound $frzDir/run_part.spec restart_file
sed -i "/restart_file /crestart_file output/${caseName}/restart.nc" $frzDir/run_part.spec
sed -i "/t_max /ct_max $3                           # total simulation time (s)" $frzDir/run_part.spec
comment_for_restart $frzDir/run_part.spec
sleep 3

../../build/partmc run_part.spec
uncomment_for_restart run_part.spec

sed -i "/restart /crestart no                      # whether to restart from saved state (yes/no)" $frzDir/run_part.spec
add_pound $frzDir/run_part.spec restart_file
