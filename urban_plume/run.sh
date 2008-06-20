#!/bin/sh

cat <<ENDINFO

Urban Plume Test-case
---------------------

This simulates an idealized urban plume with gas and aerosol
chemistry.

ENDINFO
sleep 1

#echo ../src/partmc urban_plume_with_coag.spec
#../src/partmc urban_plume_with_coag.spec
echo ../src/partmc urban_plume_no_coag.spec
../src/partmc urban_plume_no_coag.spec
