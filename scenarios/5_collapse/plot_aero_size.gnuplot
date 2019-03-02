# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set key top left

set title "Aerosol size distributions (with coag)"

set logscale

set xlabel "diameter / um"
set xrange [0.01:10]

set ytics nomirror
set y2tics

set multiplot layout 2,1

set ylabel "number concentration / (#/cm^3)"

plot "out/ogo_test_0001_aero_size_num.txt" using ($1*1e6):($2/1e6) axes x1y1 with lines title "0 hours", \
     "out/ogo_test_0001_aero_size_num.txt" using ($1*1e6):($14/1e6) axes x1y1 with lines title "12 hours", \
     "out/ogo_test_0001_aero_size_num.txt" using ($1*1e6):($32/1e6) axes x1y1 with lines title "30 hours", \
     "out/ogo_test_0001_aero_size_num.txt" using ($1*1e6):($56/1e6) axes x1y1 with lines title "54 hours", \
     "out/ogo_test_0001_aero_size_num.txt" using ($1*1e6):($74/1e6) axes x1y1 with lines title "72 hours"

set ylabel "mass concentration / (ug/m^3)"

plot "out/ogo_test_0001_aero_size_mass.txt" using ($1*1e6):($2*1e9) axes x1y1 with lines title "0 hours", \
     "out/ogo_test_0001_aero_size_mass.txt" using ($1*1e6):($14*1e9) axes x1y1 with lines title "12 hours", \
     "out/ogo_test_0001_aero_size_mass.txt" using ($1*1e6):($32*1e9) axes x1y1 with lines title "30 hours", \
     "out/ogo_test_0001_aero_size_mass.txt" using ($1*1e6):($56*1e9) axes x1y1 with lines title "54 hours", \
     "out/ogo_test_0001_aero_size_mass.txt" using ($1*1e6):($74*1e9) axes x1y1 with lines title "72 hours"

unset multiplot
