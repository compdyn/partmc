#!/bin/bash

N_ensemble=10
for expName in exp1 exp2 exp3 exp4 exp5 exp6 exp7 exp8
do
    echo "Processing $expName ..."

    for i_ensemble in $(seq 1 $N_ensemble); do
        prefix=$(printf "freezing_part_%04d " $i_ensemble)
        ../../build/freezing_process out/$expName/$prefix
    done
    
done
