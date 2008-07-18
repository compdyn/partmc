#!/bin/sh

cat <<ENDINFO

Urban Plume Test-case
---------------------

This simulates the urban plume test case for different NOx/NH3 emission szenarios.

ENDINFO
sleep 1

#echo ../src/partmc urban_plume_NOxNH3_0.5.spec
#../src/partmc urban_plume_NOxNH3_0.5.spec
#echo ../src/partmc urban_plume_NOxNH3_0.25.spec
#../src/partmc urban_plume_NOxNH3_0.25.spec
echo ../src/partmc urban_plume_no_coag.spec
../src/partmc urban_plume_no_coag.spec
#echo ../src/partmc urban_plume_with_coag.spec
#../src/partmc urban_plume_with_coag.spec

