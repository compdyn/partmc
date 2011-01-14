# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set key top left

set title "Aerosol bulk mass concentrations (no coag)"

set xrange [0:48]
set xtics 3

set xlabel "time (hours)"
set ylabel "mass concentration (ug/m^3)"

#    column  1: time (s)
#    column  2: SO4 (kg/m^3)
#    column  3: NO3 (kg/m^3)
#    column  4: Cl (kg/m^3)
#    column  5: NH4 (kg/m^3)
#    column  6: MSA (kg/m^3)
#    column  7: ARO1 (kg/m^3)
#    column  8: ARO2 (kg/m^3)
#    column  9: ALK1 (kg/m^3)
#    column 10: OLE1 (kg/m^3)
#    column 11: API1 (kg/m^3)
#    column 12: API2 (kg/m^3)
#    column 13: LIM1 (kg/m^3)
#    column 14: LIM2 (kg/m^3)
#    column 15: CO3 (kg/m^3)
#    column 16: Na (kg/m^3)
#    column 17: Ca (kg/m^3)
#    column 18: OIN (kg/m^3)
#    column 19: OC (kg/m^3)
#    column 20: BC (kg/m^3)
#    column 21: H2O (kg/m^3)

set ytics nomirror
set y2tics

set multiplot layout 2,1

plot "out/urban_plume2_nc_aero_species.txt" using ($1/3600):($3*1e9) axes x1y1 with lines title "NO3", \
     "out/urban_plume2_nc_aero_species.txt" using ($1/3600):($19*1e9) axes x1y1 with lines title "OC", \
     "out/urban_plume2_nc_aero_species.txt" using ($1/3600):($5*1e9) axes x1y1 with lines title "NH4"

plot "out/urban_plume2_nc_aero_species.txt" using ($1/3600):($2*1e9) axes x1y1 with lines title "SO4", \
     "out/urban_plume2_nc_aero_species.txt" using ($1/3600):($20*1e9) axes x1y1 with lines title "BC", \
     "out/urban_plume2_nc_aero_species.txt" using ($1/3600):(($7+$8+$9+$10)*1e9) axes x1y1 with lines title "SOA" # ARO1 + ARO2 + ALK1 + OLE1

unset multiplot
