#!/bin/sh

cat <<ENDINFO

Urban Plume Test-case
---------------------

This simulates an idealized urban plume with gas and aerosol
chemistry.

ENDINFO
sleep 1

echo ../build/partmc urban_plume_with_coag.spec
../build/partmc urban_plume_with_coag.spec
echo ../build/partmc urban_plume_no_coag.spec
../build/partmc urban_plume_no_coag.spec
