#!/bin/sh

cat <<ENDINFO

Urban Plume Test-case
---------------------

This simulates an idealized urban plume with gas and aerosol
chemistry.

ENDINFO
sleep 1

echo ../src/partmc urban_plume_test.spec
../src/partmc urban_plume_test.spec

./plot_gas.py
./plot_env.py
./plot_aero_hist.py
./plot_aero_dist_mass.py
./plot_aero_comp_bc.py
