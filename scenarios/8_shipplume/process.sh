#!/bin/sh

cd out_redist_HC_from2pm

../../../build/extract_aero_size_num 1e-9 1e-3 300 ship_plume_wc_0001_ ship_plume_wc_0001_aero_size_num.txt
gnuplot ../plot_aero_size_num.gnuplot
epstopdf aero_size_num.eps

../../../build/extract_aero_size_mass 1e-9 1e-3 300 ship_plume_wc_0001_ ship_plume_wc_0001_aero_size_mass.txt
gnuplot ../plot_aero_size_mass.gnuplot
epstopdf aero_size_mass.eps

../../../build/extract_aero_total ship_plume_wc_0001_ ship_plume_wc_0001_aero_total.txt
gnuplot ../plot_aero_total.gnuplot
epstopdf aero_total.eps

../../../build/extract_aero_species ship_plume_wc_0001_ ship_plume_wc_0001_aero_species.txt
gnuplot ../plot_aero_species.gnuplot
epstopdf aero_species.eps

../../../build/extract_env ship_plume_wc_0001_ ship_plume_wc_0001_env.txt
gnuplot ../plot_env.gnuplot
epstopdf env.eps

../../../build/extract_gas ship_plume_wc_0001_ ship_plume_wc_0001_gas.txt
gnuplot ../plot_gas.gnuplot
epstopdf gas.eps
