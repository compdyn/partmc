#!/bin/sh

cat <<ENDINFO

AIDA Chamber Test-case
---------------------

This simulates a chamber study with aerosol
coagulation.

ENDINFO
sleep 1

echo ../../build/partmc aida_self_preserving.spec
../../build/partmc aida_self_preserving.spec

