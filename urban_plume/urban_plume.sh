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
echo ../src/process_summary out/urban_plume_summary.d
../src/process_summary out/urban_plume_summary.d
echo 'Run finished, use plot_*.sh to generate plots'
