#!/bin/bash

# make sure that the current directory is the one where this script is
cd ${0%/*}

echo "../../partmc run_part.spec"
../../partmc run_part.spec

echo "../../extract_env out/condense_0001_ out/condense_env.txt"
../../extract_env out/condense_0001_ out/condense_env.txt
echo "../../numeric_diff true_env.txt out/condense_env.txt 0 1e-8"
../../numeric_diff true_env.txt out/condense_env.txt 0 1e-8
exit $?
