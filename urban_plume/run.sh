#!/bin/sh

cat <<ENDINFO

Urban Plume Test-case
---------------------

This simulates an idealized urban plume with gas and aerosol
chemistry.

ENDINFO
sleep 1

echo ../src/partmc urban_plume.spec
../src/partmc urban_plume.spec

SUBDIR=.

echo ./plot_aero_dist.py $SUBDIR
./plot_aero_dist.py $SUBDIR
echo ./plot_aero_hist.py $SUBDIR
./plot_aero_hist.py $SUBDIR
echo ./plot_gas.py $SUBDIR
./plot_gas.py $SUBDIR
echo ./plot_aero_comp_bc.py $SUBDIR
./plot_aero_comp_bc.py $SUBDIR
echo "View out/$SUBDIR"'*.pdf'
