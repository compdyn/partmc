#!/bin/sh

cat <<ENDINFO

Emissions and Background Dilution Test-case
-------------------------------------------

The initial condition is mono-disperse and only has a single
species. Emissions are also mono-disperse, but of a different
size. The background is a third mono-disperse distribution. The system
limits to a steady-state consisting only of the emissions and
background sizes.

This test-case uses both the particle and section codes with only
emissions and background dilution (no coagulation or
condensation). With only emissions and dilution the number
distribution n(r,t) satisfies:

 d n(r,t)
---------- = k_emit * n_emit(r) + k_dilute * (n_back(r) - n(r,t))
    dt

n(r,0) = n_init(r)

This is a family of ODEs parameterized by r with solution:

n(r,t) = (n_init(r) - n_lim(r)) * exp(-k_dilute * t) + n_lim(r)

where the steady state limit is:

                                    k_emit
n(r,inf) = n_lim(r) = n_back(r) + ---------- n_emit(r)
                                   k_dilute

ENDINFO
sleep 1

echo ../src/partmc emission_mc.spec
../src/partmc emission_mc.spec
echo ../src/partmc emission_sect.spec
../src/partmc emission_sect.spec
echo ../src/partmc emission_exact.spec
../src/partmc emission_exact.spec

echo ./emission_plot.py
./emission_plot.py
echo Now view out/emission_dist.pdf and out/emission_history.pdf
