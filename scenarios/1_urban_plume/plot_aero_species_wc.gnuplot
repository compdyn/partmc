# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set key top left

set title "Aerosol bulk mass concentrations (with coag)"

set xrange [0:24]
set xtics 3

set xlabel "time / h"
set ylabel "mass concentration / (ug/m^3)"

#    column  1: time (s)
#    column  2: aerosol number concentration (#/m^3)
#    column  3: aerosol mass concentration (kg/m^3)
#    column  4: aerosol SO4 concentration (kg/m^3)
#    column  5: aerosol NO3 concentration (kg/m^3)
#    column  6: aerosol Cl concentration (kg/m^3)
#    column  7: aerosol NH4 concentration (kg/m^3)
#    column  8: aerosol MSA concentration (kg/m^3)
#    column  9: aerosol ARO1 concentration (kg/m^3)
#    column 10: aerosol ARO2 concentration (kg/m^3)
#    column 11: aerosol ALK1 concentration (kg/m^3)
#    column 12: aerosol OLE1 concentration (kg/m^3)
#    column 13: aerosol API1 concentration (kg/m^3)
#    column 14: aerosol API2 concentration (kg/m^3)
#    column 15: aerosol LIM1 concentration (kg/m^3)
#    column 16: aerosol LIM2 concentration (kg/m^3)
#    column 17: aerosol CO3 concentration (kg/m^3)
#    column 18: aerosol Na concentration (kg/m^3)
#    column 19: aerosol Ca concentration (kg/m^3)
#    column 20: aerosol OIN concentration (kg/m^3)
#    column 21: aerosol OC concentration (kg/m^3)
#    column 22: aerosol BC concentration (kg/m^3)
#    column 23: aerosol H2O concentration (kg/m^3)

set ytics nomirror
set y2tics

set multiplot layout 2,1

plot "out/urban_plume_wc_0001_aero_time.txt" using ($1/3600):($5*1e9) axes x1y1 with lines title "NO3", \
     "out/urban_plume_wc_0001_aero_time.txt" using ($1/3600):($21*1e9) axes x1y1 with lines title "OC", \
     "out/urban_plume_wc_0001_aero_time.txt" using ($1/3600):($7*1e9) axes x1y1 with lines title "NH4"

plot "out/urban_plume_wc_0001_aero_time.txt" using ($1/3600):($4*1e9) axes x1y1 with lines title "SO4", \
     "out/urban_plume_wc_0001_aero_time.txt" using ($1/3600):($22*1e9) axes x1y1 with lines title "BC", \
     "out/urban_plume_wc_0001_aero_time.txt" using ($1/3600):(($9+$10+$11+$12+$13+$14+$15+$16)*1e9) axes x1y1 with lines title "SOA" # ARO1 + ARO2 + ALK1 + OLE1 + API1 + API2 + LIM1 + LIM2

unset multiplot
