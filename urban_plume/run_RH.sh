#!/bin/sh

cat <<ENDINFO

Urban Plume Test-case
---------------------

This simulates the urban plume test case for different RH.

ENDINFO
sleep 1

echo ../src/partmc urban_plume_RH5.spec
../src/partmc urban_plume_RH5.spec
echo ../src/partmc urban_plume_RH55.spec
../src/partmc urban_plume_RH55.spec
echo ../src/partmc urban_plume_RH90.spec
../src/partmc urban_plume_RH90.spec
