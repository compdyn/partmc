# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set key top right

set title "Environment data (with coag)"

set xrange [0:48]
set xtics 3

#  column 1: time (s)
#  column 2: temperature (K)
#  column 3: relative_humidity (1)
#  column 4: pressure (Pa)
#  column 5: mixing height (m)

set multiplot layout 2,1

set xlabel "time (hours)"
set ylabel "temperature (K)"
set y2label "relative humidity (%)"

set ytics nomirror
set y2tics

plot "out/urban_plume2_wc_env.txt" using ($1/3600):2 axes x1y1 with lines title "temperature", \
     "out/urban_plume2_wc_env.txt" using ($1/3600):($3*100) axes x1y2 with lines title "relative humidity"

set ylabel "mixing height (m)"
unset y2label

set ytics mirror
unset y2tics

plot "out/urban_plume2_wc_env.txt" using ($1/3600):5 axes x1y1 with lines title "mixing height"

unset multiplot
