#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

# The data should have already been processed by ./process.sh

gnuplot -persist plot_size_mass.gnuplot
gnuplot -persist plot_size_num.gnuplot
gnuplot -persist plot_time_mass.gnuplot
gnuplot -persist plot_time_num.gnuplot
