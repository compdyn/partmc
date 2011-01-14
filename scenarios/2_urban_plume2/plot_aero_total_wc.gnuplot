# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set key top left

set title "Aerosol bulk concentrations (with coag)"

set xrange [0:48]

set xlabel "time (hours)"
set ylabel "number concentration (#/m^3)"
set y2label "mass concentration (kg/m^3)"

set ytics nomirror
set y2tics

#    column 1: time (s)
#    column 2 = aerosol number concentration (#/m^3)
#    column 3 = aerosol mass concentration (kg/m^3)

plot "out/urban_plume2_wc_aero_total.txt" using ($1/3600):2 axes x1y1 with lines title "number concentration", \
     "out/urban_plume2_wc_aero_total.txt" using ($1/3600):3 axes x1y2 with lines title "mass concentration"
