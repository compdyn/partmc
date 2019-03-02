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

plot "out/collapse_0001_aero_size_num.txt" using ($1*1e6):($2/1e6) axes x1y1 with lines title "0 hours", \
     "out/collapse_0001_aero_size_num.txt" using ($1*1e6):($8/1e6) axes x1y1 with lines title "6 hours", \
     "out/collapse_0001_aero_size_num.txt" using ($1*1e6):($14/1e6) axes x1y1 with lines title "12 hours", \
     "out/collapse_0001_aero_size_num.txt" using ($1*1e6):($26/1e6) axes x1y1 with lines title "24 hours", \

set ylabel "mass concentration / (ug/m^3)"

plot "out/collapse_0001_aero_size_mass.txt" using ($1*1e6):($2*1e9) axes x1y1 with lines title "0 hours", \
     "out/collapse_0001_aero_size_mass.txt" using ($1*1e6):($8*1e9) axes x1y1 with lines title "6 hours", \
     "out/collapse_0001_aero_size_mass.txt" using ($1*1e6):($14*1e9) axes x1y1 with lines title "12 hours", \
     "out/collapse_0001_aero_size_mass.txt" using ($1*1e6):($26*1e9) axes x1y1 with lines title "24 hours", \

unset multiplot
