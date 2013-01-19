# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xrange [0.01:1]
set yrange [0:80]

set logscale x

set xlabel "diameter / um"
set ylabel "BC dry mass fraction / %"

#    column  1: particle ID number
#    column  2: number concentration (m^{-3})
#    column  3: particle diameter (m)
#    column  4: particle total mass (kg)
#    column  5: SO4 mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column  6: NO3 mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column  7: Cl mass (kg) - density = 0.2200E+04 (kg/m^3)
#    column  8: NH4 mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column  9: MSA mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 10: ARO1 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 11: ARO2 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 12: ALK1 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 13: OLE1 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 14: API1 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 15: API2 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 16: LIM1 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 17: LIM2 mass (kg) - density = 0.1400E+04 (kg/m^3)
#    column 18: CO3 mass (kg) - density = 0.2600E+04 (kg/m^3)
#    column 19: Na mass (kg) - density = 0.2200E+04 (kg/m^3)
#    column 20: Ca mass (kg) - density = 0.2600E+04 (kg/m^3)
#    column 21: OIN mass (kg) - density = 0.2600E+04 (kg/m^3)
#    column 22: OC mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 23: BC mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 24: H2O mass (kg) - density = 0.1000E+04 (kg/m^3)

set multiplot layout 2,2

set title "BC composition (with coag) at 0 hours"
plot "out/urban_plume2_wc_0001_00000001_aero_particles.txt" using ($3*1e6):(($23/($4-$24))*100) with points notitle
set title "BC composition (with coag) at 5 hours"
plot "out/urban_plume2_wc_0001_00000006_aero_particles.txt" using ($3*1e6):(($23/($4-$24))*100) with points notitle
set title "BC composition (with coag) at 7 hours"
plot "out/urban_plume2_wc_0001_00000008_aero_particles.txt" using ($3*1e6):(($23/($4-$24))*100) with points notitle
set title "BC composition (with coag) at 24 hours"
plot "out/urban_plume2_wc_0001_00000025_aero_particles.txt" using ($3*1e6):(($23/($4-$24))*100) with points notitle

unset multiplot
