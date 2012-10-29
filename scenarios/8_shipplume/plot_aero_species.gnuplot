
set terminal postscript eps color
set output "aero_species.eps"

set key top left

set xlabel "time / s"
set ylabel "mass conc  / (ug m^{-3})"
set logscale y

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

plot "ship_plume_wc_0001_aero_species.txt" using 1:($3*1e9) axes x1y1 with lines title "NO3", \
     "ship_plume_wc_0001_aero_species.txt" using 1:($19*1e9) axes x1y1 with lines title "OC", \
     "ship_plume_wc_0001_aero_species.txt" using 1:($5*1e9) axes x1y1 with lines title "NH4", \
     "ship_plume_wc_0001_aero_species.txt" using 1:($2*1e9) axes x1y1 with lines title "SO4", \
     "ship_plume_wc_0001_aero_species.txt" using 1:($20*1e9) axes x1y1 with lines title "BC", \
     "ship_plume_wc_0001_aero_species.txt" using 1:(($7+$8+$9+$10+$11+$12+$13+$14)*1e9) axes x1y1 with lines title "SOA" # ARO1 + ARO2 + ALK1 + OLE1 + API1 + API2 + LIM1 + LIM2
