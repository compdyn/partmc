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
echo ../src/process_summary out/emission_mc_summary.d
../src/process_summary out/emission_mc_summary.d
echo ../test/emission_summary_to_history out/emission_mc_summary.d
../test/emission_summary_to_history out/emission_mc_summary.d

echo ../src/partmc emission_sect.spec
../src/partmc emission_sect.spec
echo ../src/process_summary out/emission_sect_summary.d
../src/process_summary out/emission_sect_summary.d
echo ../test/emission_summary_to_history out/emission_sect_summary.d
../test/emission_summary_to_history out/emission_sect_summary.d

echo ../src/partmc emission_exact.spec
../src/partmc emission_exact.spec
echo ../src/process_summary out/emission_exact_summary.d
../src/process_summary out/emission_exact_summary.d
echo ../test/emission_summary_to_history out/emission_exact_summary.d
../test/emission_summary_to_history out/emission_exact_summary.d

echo Plotting number density
gnuplot -persist <<ENDNUM
set logscale x
set xlabel "radius (m)"
set ylabel "number density (#/m^3)"
set title "Emissions and background dilution"
plot [1e-5:1e-4] "out/emission_mc_summary_aero_binned.d" index 0 using 1:2 title "Monte Carlo (0 mins)"
replot "out/emission_mc_summary_aero_binned.d" index 5 using 1:2 title "Monte Carlo (5 mins)"
replot "out/emission_mc_summary_aero_binned.d" index 10 using 1:2 title "Monte Carlo (10 mins)"
replot "out/emission_sect_summary_aero_binned.d" index 0 using 1:2 w l title "Sectional (0 mins)"
replot "out/emission_sect_summary_aero_binned.d" index 5 using 1:2 w l title "Sectional (5 mins)"
replot "out/emission_sect_summary_aero_binned.d" index 10 using 1:2 w l title "Sectional (10 mins)"
replot "out/emission_exact_summary_aero_binned.d" index 0 using 1:2 w l title "Exact (0 mins)"
replot "out/emission_exact_summary_aero_binned.d" index 5 using 1:2 w l title "Exact (5 mins)"
replot "out/emission_exact_summary_aero_binned.d" index 10 using 1:2 w l title "Exact (10 mins)"
set terminal postscript eps
set output "out/emission_plot_num.eps"
replot
ENDNUM
epstopdf out/emission_plot_num.eps

echo Plotting number density time history
gnuplot -persist <<ENDHIST
set xlabel "time (s)"
set ylabel "number density (#/m^3)"
set title "Emissions and background dilution"
plot "out/emission_mc_summary_history.d" using 1:2 title "init cond (MC)"
replot "out/emission_mc_summary_history.d" using 1:3 title "emissions (MC)"
replot "out/emission_mc_summary_history.d" using 1:4 title "background (MC)"
replot "out/emission_sect_summary_history.d" using 1:2 w l title "init cond (sect)"
replot "out/emission_sect_summary_history.d" using 1:3 w l title "emissions (sect)"
replot "out/emission_sect_summary_history.d" using 1:4 w l title "background (sect)"
replot "out/emission_exact_summary_history.d" using 1:2 w l title "init cond (exact)"
replot "out/emission_exact_summary_history.d" using 1:3 w l title "emissions (exact)"
replot "out/emission_exact_summary_history.d" using 1:4 w l title "background (exact)"
set terminal postscript eps
set output "out/emission_plot_num_history.eps"
replot
ENDHIST
epstopdf out/emission_plot_num_history.eps
