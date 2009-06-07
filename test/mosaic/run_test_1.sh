#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

# Poisson mean m, prob of k events is:
# f(k,m) = m^k e^(-m) / k!
# so f(3,3) = e^(-1) = 0.22
# run below will fail if we don't generate exactly three particles
# so prob(failure) = 1 - f(1,1) = 0.78
# using 40 tries, prob(failure on all 40 tries) = 0.78^40 = 4e-5

mc_run_not_ok=1
MAX_TRIES=40
try_number=0
while (( $mc_run_not_ok )) ; do
    try_number=$(( $try_number + 1 ))
    echo Try number $try_number
    echo "../../partmc run_mc.spec"
    ../../partmc run_mc.spec
    echo "../../extract_state_aero_size_num 1e-8 1e-3 160 out/mosaic_0001_ out/mosaic_aero_size_num.txt"
    ../../extract_state_aero_size_num 1e-8 1e-3 160 out/mosaic_0001_ out/mosaic_aero_size_num.txt
    echo "../../numeric_diff true_aero_size_num.txt out/mosaic_aero_size_num.txt 0 1e-8 0 0 1 2"
    ../../numeric_diff true_aero_size_num.txt out/mosaic_aero_size_num.txt 0 1e-8 0 0 1 2
    mc_run_not_ok=$?
    if (( $mc_run_not_ok )) ; then
	if (( $try_number > $MAX_TRIES )) ; then
	    echo "Maximum number of tries exceeded: giving up..."
	    exit 1
	fi
	echo "Retrying..."
    fi
done

echo "../../extract_state_aero_size_num 1e-8 1e-3 160 out/mosaic_0001_ out/mosaic_aero_size_num.txt"
../../extract_state_aero_size_num 1e-8 1e-3 160 out/mosaic_0001_ out/mosaic_aero_size_num.txt
echo "../../numeric_diff true_aero_size_num.txt out/mosaic_aero_size_num.txt 0 1e-8 0 0 0 0"
../../numeric_diff true_aero_size_num.txt out/mosaic_aero_size_num.txt 0 1e-8 0 0 0 0
exit $?
