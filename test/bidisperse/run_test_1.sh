#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

# Poisson mean m, prob of k events is:
# f(k,m) = m^k e^(-m) / k!
# so f(1,1) = e^(-1) = 0.36
# run below will fail if we don't generate exactly one large particle
# so prob(failure) = 1 - f(1,1) = 0.63
# using 20 tries, prob(failure on all 20 tries) = 0.63^20 = 1e-4

mc_run_not_ok=1
MAX_TRIES=20
try_number=0
while (( $mc_run_not_ok )) ; do
    try_number=$(( $try_number + 1 ))
    echo Try number $try_number
    ../../partmc run_part.spec
    if ! ../../test_bidisperse_extract ; then
	if (( $try_number > $MAX_TRIES )) ; then
	    echo "Maximum number of tries exceeded: giving up..."
	    exit 1
	fi
	echo "Retrying..."
    else
        mc_run_not_ok=0
    fi
done

../../test_bidisperse_ode

# extract size distributions for plotting
../../extract_aero_size_num 1e-5 1e-3 255 out/bidisperse_part_0001_ out/bidisperse_part_aero_size_num.txt

../../numeric_diff out/bidisperse_part_data.txt out/bidisperse_ode_data.txt 0 1e-5 0 0 1 1
exit $?
