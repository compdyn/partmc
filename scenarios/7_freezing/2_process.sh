#!/bin/bash

N_ensemble=10

for i_ensemble in $(seq 1 $N_ensemble); do
	prefix=$(printf "freezing_part_%04d " $i_ensemble)
	../../build/freezing_process out/$prefix
done
    
