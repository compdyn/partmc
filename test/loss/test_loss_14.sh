#!/bin/bash

set -e
set -v
cd ${0%/*}
mkdir -p out

if ! grep -q '^ENABLE_QUADPACK:BOOL=ON$' ../../CMakeCache.txt ; then
	echo "Skipping QUADPACK modal drydep test in non-QUADPACK build"
	exit 0
fi

../../partmc run_modal_drydep.spec
../../extract_sectional_aero_time out/loss_modal_drydep

../../numeric_diff --by col --rel-tol 1e-8 \
		loss_modal_drydep_aero_time_ref.txt \
		out/loss_modal_drydep_aero_time.txt