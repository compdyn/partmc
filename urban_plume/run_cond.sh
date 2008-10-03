#!/bin/sh

cat <<ENDINFO

Urban Plume Test-case
---------------------

This simulates an idealized urban plume with gas and aerosol
chemistry.

ENDINFO
sleep 1

echo ../src/partmc cond_test.spec
../src/partmc cond_test.spec
