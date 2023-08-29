#!/bin/bash

expName=chiexp
t_max=3600
frzDir=/data/keeling/a/wenhant2/modeling/partmc/scenarios/7_freezing/
#mkdir -p $frzDir/output/$expName
#for chi in 0.2 0.4 0.6 0.8
for chi in 0.001 0.01 0.05 0.1 0.15
#for chi in 0.999
do
    caseName=${expName}_${chi}
    echo $caseName
    ./freezing_run_chi.sh $caseName $chi $t_max
done

