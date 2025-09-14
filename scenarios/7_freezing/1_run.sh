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

# Copy the constant settings among experiments
cp -p $setDir/*.dat .
cp -p $setDir/run_part.spec .
echo "Running $caseName ..."
# Copy the special settings for that experiment
cp -p $setDir/$expName/*.dat .
sed -i "/output_prefix /coutput_prefix ${outDir}/freezing_part # prefix of output files" run_part.spec
mkdir -p $outDir
cp -p run_part.spec $outDir
cp -p *.dat $outDir
sleep 1

../../build/partmc run_part.spec
    
