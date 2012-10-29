#!/bin/sh

cat <<ENDINFO

Ship Plume Test-case
---------------------

This simulates a ship plume with gas and aerosol
chemistry.

ENDINFO
sleep 1

echo ../../build/partmc ship_plume_no_coag.spec
../../build/partmc ship_plume_no_coag.spec
