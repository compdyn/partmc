#!/bin/sh

# exit on error
set -e
# turn on command echoing
set -v

# The data should have already been processed by ./2_process.sh

gnuplot -persist plot_size_num_linear.gnuplot
gnuplot -persist plot_size_num_log.gnuplot 
gnuplot -persist plot_size_mass_linear.gnuplot 
gnuplot -persist plot_size_mass_log.gnuplot
