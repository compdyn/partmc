# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale
set xlabel "diameter / m"
set ylabel "number concentration / (1/m^3)"

set xrange [1e-9:1e-5]
set yrange [1e1:1e10]

plot "out/weighting_flat_0001_aero_size_num.txt" using 1:2 title "flat t = 0 hours", \
     "out/weighting_flat_specified_0001_aero_size_num.txt" using 1:2 title "flat specified t = 0 hours", \
     "out/weighting_flat_source_0001_aero_size_num.txt" using 1:2 title "flat source t = 0 hours", \
     "out/weighting_exact_aero_size_num.txt" using 1:2 with lines title "sectional t = 0 hours", \
