# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "radius (m)"
set ylabel "number concentration (1/m^3)"

set xrange [1e-9:1e-6]
set yrange [0:5e8]

plot "out/average_size_num.txt" using 1:2 with lines title "base", \
     "out/average_comp_size_num.txt" using 1:2 with points pointsize 3 title "comp average", \
     "out/average_sizenum_size_num.txt" using 1:2 with points pointsize 2 title "size(num) average", \
     "out/average_sizevol_size_num.txt" using 1:2 with lines title "size(vol) average", \
     "out/average_compsizenum_size_num.txt" using 1:2 with points pointsize 1 pointtype 4 title "comp then size(num) average", \
     "out/average_compsizevol_size_num.txt" using 1:2 with points title "comp then size(vol) average"
