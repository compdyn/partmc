#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

mkdir -p out

../../build/partmc urban_plume_aq_chem.spec
../../build/partmc urban_plume_aq_chem_2.spec

# Now run ./process.sh to process the data
