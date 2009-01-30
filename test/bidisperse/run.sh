#!/bin/sh

cat README
sleep 1

echo ../../build/partmc run_mc.spec
../../build/partmc run_mc.spec
echo ../../build/bidisperse_ode
../../build/bidisperse_ode

echo ./plot.py
./plot.py
echo "Now view out/bidisperse_*.pdf"
