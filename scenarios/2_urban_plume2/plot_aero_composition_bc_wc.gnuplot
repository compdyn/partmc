# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xrange [0.01:1]
set yrange [0:80]

set logscale x

set xlabel "diameter (um)"
set ylabel "BC dry mass fraction (%)"

#    column  1: particle ID number
#    column  2: computational volume (m^3)
#    column  3: particle radius (m)
#    column  4: particle total mass (kg)
#    column  5: SO4 mass (kg) - density =  0.180E+04 (kg/m^3)
#    column  6: NO3 mass (kg) - density =  0.180E+04 (kg/m^3)
#    column  7: Cl mass (kg) - density =  0.220E+04 (kg/m^3)
#    column  8: NH4 mass (kg) - density =  0.180E+04 (kg/m^3)
#    column  9: MSA mass (kg) - density =  0.180E+04 (kg/m^3)
#    column 10: ARO1 mass (kg) - density =  0.140E+04 (kg/m^3)
#    column 11: ARO2 mass (kg) - density =  0.140E+04 (kg/m^3)
#    column 12: ALK1 mass (kg) - density =  0.140E+04 (kg/m^3)
#    column 13: OLE1 mass (kg) - density =  0.140E+04 (kg/m^3)
#    column 14: API1 mass (kg) - density =  0.140E+04 (kg/m^3)
#    column 15: API2 mass (kg) - density =  0.140E+04 (kg/m^3)
#    column 16: LIM1 mass (kg) - density =  0.140E+04 (kg/m^3)
#    column 17: LIM2 mass (kg) - density =  0.140E+04 (kg/m^3)
#    column 18: CO3 mass (kg) - density =  0.260E+04 (kg/m^3)
#    column 19: Na mass (kg) - density =  0.220E+04 (kg/m^3)
#    column 20: Ca mass (kg) - density =  0.260E+04 (kg/m^3)
#    column 21: OIN mass (kg) - density =  0.260E+04 (kg/m^3)
#    column 22: OC mass (kg) - density =  0.140E+04 (kg/m^3)
#    column 23: BC mass (kg) - density =  0.180E+04 (kg/m^3)
#    column 24: H2O mass (kg) - density =  0.100E+04 (kg/m^3)

set multiplot layout 2,2

set title "BC composition (with coag) at 0 hours"
plot "out/urban_plume_wc_aero_particle_mass_00000001.txt" using ($3*2e6):(($23/($4-$24))*100) with points notitle
set title "BC composition (with coag) at 5 hours"
plot "out/urban_plume_wc_aero_particle_mass_00000006.txt" using ($3*2e6):(($23/($4-$24))*100) with points notitle
set title "BC composition (with coag) at 7 hours"
plot "out/urban_plume_wc_aero_particle_mass_00000008.txt" using ($3*2e6):(($23/($4-$24))*100) with points notitle
set title "BC composition (with coag) at 24 hours"
plot "out/urban_plume_wc_aero_particle_mass_00000025.txt" using ($3*2e6):(($23/($4-$24))*100) with points notitle

unset multiplot
