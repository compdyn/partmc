#!/bin/sh

cat README
sleep 1

echo ../../src/partmc run_mc.spec
../../src/partmc run_mc.spec
echo ./bidisperse_ode
./bidisperse_ode

echo ./plot.py
./plot.py
echo "Now view out/bidisperse_*.pdf"
