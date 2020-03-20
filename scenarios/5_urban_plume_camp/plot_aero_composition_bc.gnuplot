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
#    column  5: dust.LHD_DUST mass (kg) - density = 0.2650E+04 (kg/m^3)
#    column  6: dust.LLD_DUST mass (kg) - density = 0.2650E+04 (kg/m^3)
#    column  7: sea_salt.SEA_SALT mass (kg) - density = 0.2160E+04 (kg/m^3)
#    column  8: organic_matter.PM_phob mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column  9: organic_matter.PM_phil mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 10: organic_matter.ISOP-P1_aero mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 11: organic_matter.ISOP-P2_aero mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 12: organic_matter.TERP-P1_aero mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 13: organic_matter.TERP-P2_aero mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 14: black_carbon.BC_phob mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 15: black_carbon.BC_phil mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 16: other_PM.other_PM mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 17: other_PM.other_other_PM mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 18: aqueous.H2O_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 19: aqueous.H_p mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 20: aqueous.OH_m mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 21: aqueous.H2SO4_aq mass (kg) - density = 0.1840E+04 (kg/m^3)
#    column 22: aqueous.SO4_mm mass (kg) - density = 0.1840E+04 (kg/m^3)
#    column 23: aqueous.HSO4_m mass (kg) - density = 0.1840E+04 (kg/m^3)
#    column 24: aqueous.HNO3_aq mass (kg) - density = 0.1510E+04 (kg/m^3)
#    column 25: aqueous.NO3_m mass (kg) - density = 0.1510E+04 (kg/m^3)
#    column 26: aqueous.NH3_aq mass (kg) - density = 0.8800E+03 (kg/m^3)
#    column 27: aqueous.NH4_p mass (kg) - density = 0.8800E+03 (kg/m^3)
#    column 28: aqueous.HCL_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 29: aqueous.CL_m mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 30: aqueous.NO2_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 31: aqueous.NO_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 32: aqueous.O3_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 33: aqueous.NO3_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 34: aqueous.N2O5_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 35: aqueous.HONO_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 36: aqueous.HNO4_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 37: aqueous.H2O2_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 38: aqueous.NTR_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 39: aqueous.ROOH_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 40: aqueous.FORM_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 41: aqueous.ALD2_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 42: aqueous.ALDX_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 43: aqueous.CO_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 44: aqueous.MEPX_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 45: aqueous.MEOH_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 46: aqueous.FACD_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 47: aqueous.PAN_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 48: aqueous.PACD_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 49: aqueous.AACD_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 50: aqueous.PANX_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 51: aqueous.SO2_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 52: aqueous.CL2_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 53: aqueous.HOCL_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 54: aqueous.FMCL_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 55: aqueous.ISOP-P1_aero mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 56: aqueous.ISOP-P2_aero mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 57: aqueous.TERP-P1_aero mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 58: aqueous.TERP-P2_aero mass (kg) - density = 0.1800E+04 (kg/m^3)
#    column 59: aqueous.DMS_aq mass (kg) - density = 0.1000E+04 (kg/m^3)
#    column 60: aqueous.ETOH_aq mass (kg) - density = 0.1000E+04 (kg/m^3)

set multiplot layout 2,2

set title "BC composition (with coag) at 0 hours"
plot "out/urban_plume_0001_00000001_aero_particles.txt" using ($3*1e6):((($14+$15)/($4-$18))*100) with points notitle
set title "BC composition (with coag) at 5 hours"
plot "out/urban_plume_0001_00000006_aero_particles.txt" using ($3*1e6):((($14+$15)/($4-$18))*100) with points notitle
set title "BC composition (with coag) at 7 hours"
plot "out/urban_plume_0001_00000008_aero_particles.txt" using ($3*1e6):((($14+$15)/($4-$18))*100) with points notitle
set title "BC composition (with coag) at 24 hours"
plot "out/urban_plume_0001_00000025_aero_particles.txt" using ($3*1e6):((($14+$15)/($4-$18))*100) with points notitle

unset multiplot
