# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set key top left

#set title "Aerosol bulk concentrations (with coag)"

#set xrange [0:24]

set xlabel "time (hours)"
set ylabel "number concentration (#/m^3)"
#set y2label "mass concentration (kg/m^3)"

#set ytics nomirror
#set y2tics

#    column 1: time (s)
#    column 2 = aerosol number concentration (#/m^3)
#    column 3 = aerosol mass concentration (kg/m^3)

plot "out/parallel_0001_aero_total.txt" using 1:2 axes x1y1 with lines title "num p1", \
     "out/parallel_0002_aero_total.txt" using 1:2 axes x1y1 with lines title "num p2", \
     "out/serial_aero_total.txt" using 1:2 axes x1y1 with lines title "num serial", \
     "out/sect_aero_total.txt" using 1:2 axes x1y1 with lines title "num sect"
