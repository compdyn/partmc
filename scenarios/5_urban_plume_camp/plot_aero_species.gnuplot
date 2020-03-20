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
#    column  4: aerosol dust.LHD_DUST concentration (kg/m^3)
#    column  5: aerosol dust.LLD_DUST concentration (kg/m^3)
#    column  6: aerosol sea_salt.SEA_SALT concentration (kg/m^3)
#    column  7: aerosol organic_matter.PM_phob concentration (kg/m^3)
#    column  8: aerosol organic_matter.PM_phil concentration (kg/m^3)
#    column  9: aerosol organic_matter.ISOP-P1_aero concentration (kg/m^3)
#    column 10: aerosol organic_matter.ISOP-P2_aero concentration (kg/m^3)
#    column 11: aerosol organic_matter.TERP-P1_aero concentration (kg/m^3)
#    column 12: aerosol organic_matter.TERP-P2_aero concentration (kg/m^3)
#    column 13: aerosol black_carbon.BC_phob concentration (kg/m^3)
#    column 14: aerosol black_carbon.BC_phil concentration (kg/m^3)
#    column 15: aerosol other_PM.other_PM concentration (kg/m^3)
#    column 16: aerosol other_PM.other_other_PM concentration (kg/m^3)
#    column 17: aerosol aqueous.H2O_aq concentration (kg/m^3)
#    column 18: aerosol aqueous.H_p concentration (kg/m^3)
#    column 19: aerosol aqueous.OH_m concentration (kg/m^3)
#    column 20: aerosol aqueous.H2SO4_aq concentration (kg/m^3)
#    column 21: aerosol aqueous.SO4_mm concentration (kg/m^3)
#    column 22: aerosol aqueous.HSO4_m concentration (kg/m^3)
#    column 23: aerosol aqueous.HNO3_aq concentration (kg/m^3)
#    column 24: aerosol aqueous.NO3_m concentration (kg/m^3)
#    column 25: aerosol aqueous.NH3_aq concentration (kg/m^3)
#    column 26: aerosol aqueous.NH4_p concentration (kg/m^3)
#    column 27: aerosol aqueous.HCL_aq concentration (kg/m^3)
#    column 28: aerosol aqueous.CL_m concentration (kg/m^3)
#    column 29: aerosol aqueous.NO2_aq concentration (kg/m^3)
#    column 30: aerosol aqueous.NO_aq concentration (kg/m^3)
#    column 31: aerosol aqueous.O3_aq concentration (kg/m^3)
#    column 32: aerosol aqueous.NO3_aq concentration (kg/m^3)
#    column 33: aerosol aqueous.N2O5_aq concentration (kg/m^3)
#    column 34: aerosol aqueous.HONO_aq concentration (kg/m^3)
#    column 35: aerosol aqueous.HNO4_aq concentration (kg/m^3)
#    column 36: aerosol aqueous.H2O2_aq concentration (kg/m^3)
#    column 37: aerosol aqueous.NTR_aq concentration (kg/m^3)
#    column 38: aerosol aqueous.ROOH_aq concentration (kg/m^3)
#    column 39: aerosol aqueous.FORM_aq concentration (kg/m^3)
#    column 40: aerosol aqueous.ALD2_aq concentration (kg/m^3)
#    column 41: aerosol aqueous.ALDX_aq concentration (kg/m^3)
#    column 42: aerosol aqueous.CO_aq concentration (kg/m^3)
#    column 43: aerosol aqueous.MEPX_aq concentration (kg/m^3)
#    column 44: aerosol aqueous.MEOH_aq concentration (kg/m^3)
#    column 45: aerosol aqueous.FACD_aq concentration (kg/m^3)
#    column 46: aerosol aqueous.PAN_aq concentration (kg/m^3)
#    column 47: aerosol aqueous.PACD_aq concentration (kg/m^3)
#    column 48: aerosol aqueous.AACD_aq concentration (kg/m^3)
#    column 49: aerosol aqueous.PANX_aq concentration (kg/m^3)
#    column 50: aerosol aqueous.SO2_aq concentration (kg/m^3)
#    column 51: aerosol aqueous.CL2_aq concentration (kg/m^3)
#    column 52: aerosol aqueous.HOCL_aq concentration (kg/m^3)
#    column 53: aerosol aqueous.FMCL_aq concentration (kg/m^3)
#    column 54: aerosol aqueous.ISOP-P1_aero concentration (kg/m^3)
#    column 55: aerosol aqueous.ISOP-P2_aero concentration (kg/m^3)
#    column 56: aerosol aqueous.TERP-P1_aero concentration (kg/m^3)
#    column 57: aerosol aqueous.TERP-P2_aero concentration (kg/m^3)
#    column 58: aerosol aqueous.DMS_aq concentration (kg/m^3)
#    column 59: aerosol aqueous.ETOH_aq concentration (kg/m^3)

set ytics nomirror
set y2tics

set multiplot layout 2,1

plot "out/urban_plume_0001_aero_time.txt" using ($1/3600):(($23+$24)*1e9) axes x1y1 with lines title "NO3", \
     "out/urban_plume_0001_aero_time.txt" using ($1/3600):(($7+$8)*1e9) axes x1y1 with lines title "POA", \
     "out/urban_plume_0001_aero_time.txt" using ($1/3600):(($25+$26)*1e9) axes x1y1 with lines title "NH4", \
     "out/urban_plume_0001_aero_time.txt" using ($1/3600):(($17)*1e9) axes x1y1 with lines title "H2O"

plot "out/urban_plume_0001_aero_time.txt" using ($1/3600):(($9+$10+$11+$12+$54+$55+$56+$57)*1e9) axes x1y1 with lines title "SOA" # ISOP-P1 + ISOP-P2 + TERP-P1 + TERP-P2

unset multiplot
