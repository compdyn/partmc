#!/bin/sh

cat <<ENDINFO

Chamber case based on data from Naumann [2003]
---------------------

This simulates a chamber study with aerosol
coagulation.

ENDINFO
sleep 1

echo ../../build/partmc naumann_with_coag_df_2.spec
../../build/partmc naumann_with_coag_df_2.spec

