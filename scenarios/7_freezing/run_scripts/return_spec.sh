#!/bin/bash
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

frzDir=/data/nriemer/a/wenhant2/modeling/partmc/scenarios/7_freezing
cd $frzDir
uncomment_for_restart run_part.spec

sed -i "/restart /crestart no                      # whether to restart from saved state (yes/no)" $frzDir/run_part.spec
add_pound $frzDir/run_part.spec restart_file
