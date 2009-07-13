# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "radius (m)"
set ylabel "number concentration (#/m^3)"

set xrange [1e-9:1e-6]
#set yrange [1e7:1e11]

set terminal aqua 0
plot "out/parallel_0001_aero_size_num.txt" using 1:2 with linespoints title "p1 t = 0 hours", \
     "out/parallel_0002_aero_size_num.txt" using 1:2 with linespoints title "p2 t = 0 hours", \
     "out/serial_aero_size_num.txt" using 1:2 with linespoints title "serial t = 0 hours", \
     "out/sect_aero_size_num.txt" using 1:2 with lines title "sect t = 0 hours"

set terminal aqua 1
plot "out/parallel_0001_aero_size_num.txt" using 1:14 with linespoints title "p1 t = 12 hours", \
     "out/parallel_0002_aero_size_num.txt" using 1:14 with linespoints title "p2 t = 12 hours", \
     "out/serial_aero_size_num.txt" using 1:14 with linespoints title "serial t = 12 hours", \
     "out/sect_aero_size_num.txt" using 1:14 with lines title "sect t = 12 hours"

set terminal aqua 2
plot "out/parallel_0001_aero_size_num.txt" using 1:26 with linespoints title "p1 t = 24 hours", \
     "out/parallel_0002_aero_size_num.txt" using 1:26 with linespoints title "p2 t = 24 hours", \
     "out/serial_aero_size_num.txt" using 1:26 with linespoints title "serial t = 24 hours", \
     "out/sect_aero_size_num.txt" using 1:26 with lines title "sect t = 24 hours"
