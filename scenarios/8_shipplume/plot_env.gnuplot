
set terminal postscript eps color
set output "env.eps"

set key top right

set xlabel "time / s"
set ylabel "temperature / K"
set y2label "relative humidity / %"

set ytics nomirror
set y2tics

plot "ship_plume_wc_0001_env.txt" using ($1):($2) axes x1y1 with lines title "temperature", \
     "ship_plume_wc_0001_env.txt" using ($1):($3*100) axes x1y2 with lines title "relative humidity"
