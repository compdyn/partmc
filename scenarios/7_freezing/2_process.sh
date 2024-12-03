#!/bin/bash

for expName in exp1 exp2 exp3 exp4 exp5 exp6 exp7 exp8
do
    echo "Processing $expName ..."
    ../../build/freezing_process $expName
    
    
done
