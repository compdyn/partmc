#!/bin/bash

set -e
set -v
cd ${0%/*}
mkdir -p out

../../partmc run_drydep_sect.spec
../../extract_sectional_aero_time out/loss_drydep_sect

../../numeric_diff --by col --rel-tol 1e-8 \
		loss_drydep_sect_aero_time_ref.txt \
		out/loss_drydep_sect_aero_time.txt
