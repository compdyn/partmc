# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set grid
set xlabel "dry diameter D / um"
set ylabel "mass conc ug / m^{3}"

plot "out/barrel_wc_0001_aero_size_mass.txt" using ($1*1e6):($2/1e6*log(10)) title "t = 0 s, PartMC", \
     "out/barrel_wc_0001_aero_size_mass.txt" using ($1*1e6):($11/1e6*log(10)) title "t = 70 min, PartMC", \
     "out/barrel_wc_0001_aero_size_mass.txt" using ($1*1e6):($21/1e6*log(10)) title "t = 140 min, PartMC", \
     "ref_aero_size_mass_regrid.txt" using ($1*1e6):($2/1e6) with lines title "t = 0 s, Measurement", \
     "ref_aero_size_mass_regrid.txt" using ($1*1e6):($11/1e6) with lines title "t = 70 min, Measurement", \
     "ref_aero_size_mass_regrid.txt" using ($1*1e6):($21/1e6) with lines title "t = 140 min, Measurement"
