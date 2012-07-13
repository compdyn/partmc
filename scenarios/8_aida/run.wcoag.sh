#!/bin/sh

cat <<ENDINFO

AIDA Chamber case
---------------------

This simulates a chamber study with aerosol
coagulation.

ENDINFO
sleep 1

echo ../../build/partmc aida_with_coag.spec
../../build/partmc aida_with_coag.spec

