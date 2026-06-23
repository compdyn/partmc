#!/bin/bash

# Sectional condensational-growth test (independent of test_tchem_1..4).
#
# A small (100 nm) POA seed is exposed to a reservoir of condensable SIMPOL
# semivolatiles. As TChem condenses mass onto the bins, the sectional
# bin-remap (aero_binned_redistribute) must move the grown particles up into
# higher bins while conserving total number. This test asserts both:
#   (1) the number-mean particle grows into a higher bin, and
#   (2) total number concentration is conserved by the remap.

# exit on error
set -e
# turn on command echoing
set -v
# make sure that the current directory is the one where this script is
cd ${0%/*}
# make the output directory if it doesn't exist
mkdir -p out

../../partmc run_sect_grow_with_SIMPOL.spec
../../extract_sectional_aero_size --num out/tchem_grow_sect

# Column 1 is the bin diameter; column 2 is t=0 and the last column is t=t_max.
# The populated (peak) bin must move to a higher index, and the total number
# (sum over bins, the constant log-width cancels in the ratio) must be conserved.
awk '
{ diam[NR] = $1; n_init[NR] = $2; n_final[NR] = $NF }
END {
    for (i = 1; i <= NR; i++) {
        if (n_init[i]  > pk_init)  { pk_init  = n_init[i];  bin_init  = i }
        if (n_final[i] > pk_final) { pk_final = n_final[i]; bin_final = i }
        sum_init  += n_init[i]
        sum_final += n_final[i]
    }
    ratio = sum_final / sum_init
    printf "initial peak bin = %d (D = %.3e m)\n", bin_init,  diam[bin_init]
    printf "final   peak bin = %d (D = %.3e m)\n", bin_final, diam[bin_final]
    printf "number ratio (final/initial) = %.8f\n", ratio
    fail = 0
    if (bin_final <= bin_init) {
        print "FAIL: particles did not grow into a higher bin"
        fail = 1
    }
    if (ratio > 1.0 + 1e-6 || ratio < 1.0 - 1e-6) {
        print "FAIL: total number not conserved by the bin remap"
        fail = 1
    }
    if (fail) { exit 1 }
    print "PASS: particles grew into higher bins with number conserved"
}
' out/tchem_grow_sect_aero_size_num.txt
