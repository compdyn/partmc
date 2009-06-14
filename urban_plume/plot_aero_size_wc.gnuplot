# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot <filename>.gnuplot

set key top left

set title "Aerosol size distributions (with coag)"

set logscale

set xlabel "diameter (um)"
set xrange [0.01:1]

set ytics nomirror
set y2tics

set multiplot layout 2,1

set ylabel "number concentration (#/cm^3)"

plot "out/urban_plume_wc_aero_size_num.txt" using ($1*2e6):($2/1e6) axes x1y1 with lines title "0 hours", \
     "out/urban_plume_wc_aero_size_num.txt" using ($1*2e6):($302/1e6) axes x1y1 with lines title "5 hours", \
     "out/urban_plume_wc_aero_size_num.txt" using ($1*2e6):($422/1e6) axes x1y1 with lines title "7 hours", \
     "out/urban_plume_wc_aero_size_num.txt" using ($1*2e6):($1442/1e6) axes x1y1 with lines title "24 hours"

set ylabel "mass concentration (ug/m^3)"

plot "out/urban_plume_wc_aero_size_mass.txt" using ($1*2e6):($2*1e9) axes x1y1 with lines title "0 hours", \
     "out/urban_plume_wc_aero_size_mass.txt" using ($1*2e6):($302*1e9) axes x1y1 with lines title "5 hours", \
     "out/urban_plume_wc_aero_size_mass.txt" using ($1*2e6):($422*1e9) axes x1y1 with lines title "7 hours", \
     "out/urban_plume_wc_aero_size_mass.txt" using ($1*2e6):($1442*1e9) axes x1y1 with lines title "24 hours"

unset multiplot
