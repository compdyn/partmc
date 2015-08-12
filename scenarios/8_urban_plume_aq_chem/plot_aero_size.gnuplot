# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set terminal postscript

set key top left

set title "Aerosol size distributions"

set logscale

set xlabel "diameter / um"
set xrange [0.01:10]

set ytics nomirror
set y2tics

set multiplot layout 2,1

set ylabel "number concentration / (#/cm^3)"

plot "out/urban_plume_aq_chem_b_0001_aero_size_num.txt" using ($1*1e6):($5/1e6) axes x1y1 with lines title "20 min"
#     "out/urban_plume_aq_chem_b_0001_aero_size_num.txt" using ($1*1e6):($8/1e6) axes x1y1 with lines title "0.5 hours", \
#     "out/urban_plume_aq_chem_b_0001_aero_size_num.txt" using ($1*1e6):($13/1e6) axes x1y1 with lines title "1.0 hours", \
#     "out/urban_plume_aq_chem_b_0001_aero_size_num.txt" using ($1*1e6):($23/1e6) axes x1y1 with lines title "2.0 hours"

set ylabel "mass concentration / (ug/m^3)"

plot "out/urban_plume_aq_chem_b_0001_aero_size_mass.txt" using ($1*1e6):($5*1e9) axes x1y1 with lines title "20 min"
#     "out/urban_plume_aq_chem_b_0001_aero_size_mass.txt" using ($1*1e6):($8*1e9) axes x1y1 with lines title "0.5 hours", \
#     "out/urban_plume_aq_chem_b_0001_aero_size_mass.txt" using ($1*1e6):($13*1e9) axes x1y1 with lines title "1.0 hours", \
#     "out/urban_plume_aq_chem_b_0001_aero_size_mass.txt" using ($1*1e6):($23*1e9) axes x1y1 with lines title "2.0 hours"

unset multiplot
