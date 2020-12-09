# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set key top left

set title "Gas mixing ratios (with coag)"

set xrange [0:24]
set xtics 3

set multiplot layout 3,1

set xlabel "time / h"
set ylabel "gas mixing ratio / ppb"

#    column  1: time (s)
#    column  2: gas NO2 mixing ratio (ppb)
#    column  3: gas NO mixing ratio (ppb)
#    column  4: gas O mixing ratio (ppb)
#    column  5: gas O3 mixing ratio (ppb)
#    column  6: gas NO3 mixing ratio (ppb)
#    column  7: gas O1D mixing ratio (ppb)
#    column  8: gas OH mixing ratio (ppb)
#    column  9: gas HO2 mixing ratio (ppb)
#    column 10: gas N2O5 mixing ratio (ppb)
#    column 11: gas HNO3 mixing ratio (ppb)
#    column 12: gas HONO mixing ratio (ppb)
#    column 13: gas PNA mixing ratio (ppb)
#    column 14: gas H2O2 mixing ratio (ppb)
#    column 15: gas XO2 mixing ratio (ppb)
#    column 16: gas XO2N mixing ratio (ppb)
#    column 17: gas NTR mixing ratio (ppb)
#    column 18: gas ROOH mixing ratio (ppb)
#    column 19: gas FORM mixing ratio (ppb)
#    column 20: gas ALD2 mixing ratio (ppb)
#    column 21: gas ALDX mixing ratio (ppb)
#    column 22: gas PAR mixing ratio (ppb)
#    column 23: gas CO mixing ratio (ppb)
#    column 24: gas MEO2 mixing ratio (ppb)
#    column 25: gas MEPX mixing ratio (ppb)
#    column 26: gas MEOH mixing ratio (ppb)
#    column 27: gas HCO3 mixing ratio (ppb)
#    column 28: gas FACD mixing ratio (ppb)
#    column 29: gas C2O3 mixing ratio (ppb)
#    column 30: gas PAN mixing ratio (ppb)
#    column 31: gas PACD mixing ratio (ppb)
#    column 32: gas AACD mixing ratio (ppb)
#    column 33: gas CXO3 mixing ratio (ppb)
#    column 34: gas PANX mixing ratio (ppb)
#    column 35: gas ROR mixing ratio (ppb)
#    column 36: gas OLE mixing ratio (ppb)
#    column 37: gas ETH mixing ratio (ppb)
#    column 38: gas IOLE mixing ratio (ppb)
#    column 39: gas TOL mixing ratio (ppb)
#    column 40: gas CRES mixing ratio (ppb)
#    column 41: gas TO2 mixing ratio (ppb)
#    column 42: gas TOLRO2 mixing ratio (ppb)
#    column 43: gas OPEN mixing ratio (ppb)
#    column 44: gas CRO mixing ratio (ppb)
#    column 45: gas MGLY mixing ratio (ppb)
#    column 46: gas XYL mixing ratio (ppb)
#    column 47: gas XYLRO2 mixing ratio (ppb)
#    column 48: gas ISOP mixing ratio (ppb)
#    column 49: gas ISPD mixing ratio (ppb)
#    column 50: gas ISOPRXN mixing ratio (ppb)
#    column 51: gas TERP mixing ratio (ppb)
#    column 52: gas TRPRXN mixing ratio (ppb)
#    column 53: gas SO2 mixing ratio (ppb)
#    column 54: gas SULF mixing ratio (ppb)
#    column 55: gas SULRXN mixing ratio (ppb)
#    column 56: gas ETOH mixing ratio (ppb)
#    column 57: gas ETHA mixing ratio (ppb)
#    column 58: gas CL2 mixing ratio (ppb)
#    column 59: gas CL mixing ratio (ppb)
#    column 60: gas HOCL mixing ratio (ppb)
#    column 61: gas CLO mixing ratio (ppb)
#    column 62: gas FMCL mixing ratio (ppb)
#    column 63: gas HCL mixing ratio (ppb)
#    column 64: gas TOLNRXN mixing ratio (ppb)
#    column 65: gas TOLHRXN mixing ratio (ppb)
#    column 66: gas XYLNRXN mixing ratio (ppb)
#    column 67: gas XYLHRXN mixing ratio (ppb)
#    column 68: gas BENZENE mixing ratio (ppb)
#    column 69: gas BENZRO2 mixing ratio (ppb)
#    column 70: gas BNZNRXN mixing ratio (ppb)
#    column 71: gas BNZHRXN mixing ratio (ppb)
#    column 72: gas SESQ mixing ratio (ppb)
#    column 73: gas SESQRXN mixing ratio (ppb)
#    column 74: gas M mixing ratio (ppb)
#    column 75: gas O2 mixing ratio (ppb)
#    column 76: gas N2 mixing ratio (ppb)
#    column 77: gas H2O mixing ratio (ppb)
#    column 78: gas CH4 mixing ratio (ppb)
#    column 79: gas H2 mixing ratio (ppb)
#    column 80: gas N2O mixing ratio (ppb)
#    column 81: gas DUMMY mixing ratio (ppb)
#    column 82: gas DMS mixing ratio (ppb)
#    column 83: gas NH3 mixing ratio (ppb)
#    column 84: gas ISOP-P1 mixing ratio (ppb)
#    column 85: gas ISOP-P2 mixing ratio (ppb)
#    column 86: gas TERP-P1 mixing ratio (ppb)
#    column 87: gas TERP-P2 mixing ratio (ppb)

plot "out/urban_plume_0001_gas.txt" using ($1/3600):5 axes x1y1 with lines title "O3", \
     "out/urban_plume_0001_gas.txt" using ($1/3600):2 axes x1y1 with lines title "NO2", \
     "out/urban_plume_0001_gas.txt" using ($1/3600):3 axes x1y1 with lines title "NO"

plot "out/urban_plume_0001_gas.txt" using ($1/3600):11 axes x1y1 with lines title "HNO3", \
     "out/urban_plume_0001_gas.txt" using ($1/3600):19 axes x1y1 with lines title "HCHO", \
     "out/urban_plume_0001_gas.txt" using ($1/3600):53 axes x1y1 with lines title "SO2", \
     "out/urban_plume_0001_gas.txt" using ($1/3600):51 axes x1y1 with lines title "TERP", \
     "out/urban_plume_0001_gas.txt" using ($1/3600):49 axes x1y1 with lines title "ISOP", \
     "out/urban_plume_0001_gas.txt" using ($1/3600):83 axes x1y1 with lines title "NH3"

plot "out/urban_plume_0001_gas.txt" using ($77/3600):11 axes x1y1 with lines title "H2O"
unset multiplot
