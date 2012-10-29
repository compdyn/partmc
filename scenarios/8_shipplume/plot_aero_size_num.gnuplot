
set terminal postscript eps color
set output "aero_size_num.eps"

set logscale
set grid
set xlabel "wet diameter D / um"
set ylabel "num conc n / cm^{-3}"

plot "ship_plume_wc_0001_aero_size_num.txt" using ($1*1e6):($2/1e6) title "t = 0 s", \
     "ship_plume_wc_0001_aero_size_num.txt" using ($1*1e6):($3/1e6) title "t = 10 s", \
     "ship_plume_wc_0001_aero_size_num.txt" using ($1*1e6):($6/1e6) title "t = 40 s", \
     "ship_plume_wc_0001_aero_size_num.txt" using ($1*1e6):($11/1e6) title "t = 100 s"
