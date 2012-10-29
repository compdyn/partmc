
set terminal postscript eps color
set output "aero_total.eps"

set logscale y
set logscale y2
set grid
set xlabel "time / s"
set ylabel "num conc N / cm^{-3}"
set y2label "mass conc M / (ug m^{-3})"

set ytics nomirror
set y2tics

plot "ship_plume_wc_0001_aero_total.txt" using ($1):($2/1e6) axes x1y1 with lines title "num conc", \
     "ship_plume_wc_0001_aero_total.txt" using ($1):($3*1e9) axes x1y2 with lines title "mass conc"
