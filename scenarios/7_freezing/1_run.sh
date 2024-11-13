#!/bin/bash

setDir=standard_setting
outDir=out


if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
    echo "Create the folder $outDir"
else
    echo "The folder $outDir exists"
fi

cp -rp $setDir $outDir/


#if [ -e "output" ]; then
#    unlink output
#fi
#ln -sf $outDir output

for expName in exp1 exp2 exp3 exp4 exp5 exp6 exp7 exp8
do
    caseName=${expName}
    echo "Running $caseName ..."
    cp -p $setDir/$expName/*.dat .
    cp -p $setDir/run_part.spec .
    sed -i "/output_prefix /coutput_prefix ${outDir}/${caseName}/freezing_part # prefix of output files" run_part.spec
    mkdir -p $outDir/$caseName
    cp -p run_part.spec $outDir/$caseName
    cp -p *.dat $outDir/$caseName
    sleep 1

    ../../build/partmc run_part.spec
    #mv freezing_timing.txt output/$caseName
    
done
