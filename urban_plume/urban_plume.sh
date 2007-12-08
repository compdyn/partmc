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

echo ./plot_aero_dist.py
./plot_aero_dist.py
echo ./plot_aero_hist.py
./plot_aero_hist.py
echo ./plot_gas.py
./plot_gas.py
echo ./plot_aero_kappa.py
./plot_aero_kappa.py
echo ./plot_aero_comp_bc.py
./plot_aero_comp_bc.py
echo 'View out/*.pdf'
