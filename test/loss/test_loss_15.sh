#!/bin/bash

set -e
set -v
cd ${0%/*}
mkdir -p out

if grep -q '^ENABLE_QUADPACK:BOOL=ON$' ../../CMakeCache.txt ; then
    echo "Skipping non-QUADPACK modal process test in QUADPACK build"
    exit 0
fi

../../partmc run_modal_drydep.spec
../../drydep_modal_process out/loss_modal_drydep
