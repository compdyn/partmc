# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set key top left

set title "Gas mixing ratios (w hydr.)"

set xrange [0:24]
set xtics 3

set multiplot layout 2,1

set xlabel "time / h"
set ylabel "gas mixing ratio / ppb"

#    column  1: time (s)
#    column  2: gas H2SO4 mixing ratio (ppb)
#    column  3: gas HNO3 mixing ratio (ppb)
#    column  4: gas HCl mixing ratio (ppb)
#    column  5: gas NH3 mixing ratio (ppb)
#    column  6: gas NO mixing ratio (ppb)
#    column  7: gas NO2 mixing ratio (ppb)
#    column  8: gas NO3 mixing ratio (ppb)
#    column  9: gas N2O5 mixing ratio (ppb)
#    column 10: gas HONO mixing ratio (ppb)
#    column 11: gas HNO4 mixing ratio (ppb)
#    column 12: gas O3 mixing ratio (ppb)
#    column 13: gas O1D mixing ratio (ppb)
#    column 14: gas O3P mixing ratio (ppb)
#    column 15: gas OH mixing ratio (ppb)
#    column 16: gas HO2 mixing ratio (ppb)
#    column 17: gas H2O2 mixing ratio (ppb)
#    column 18: gas CO mixing ratio (ppb)
#    column 19: gas SO2 mixing ratio (ppb)
#    column 20: gas CH4 mixing ratio (ppb)
#    column 21: gas C2H6 mixing ratio (ppb)
#    column 22: gas CH3O2 mixing ratio (ppb)
#    column 23: gas ETHP mixing ratio (ppb)
#    column 24: gas HCHO mixing ratio (ppb)
#    column 25: gas CH3OH mixing ratio (ppb)
#    column 26: gas ANOL mixing ratio (ppb)
#    column 27: gas CH3OOH mixing ratio (ppb)
#    column 28: gas ETHOOH mixing ratio (ppb)
#    column 29: gas ALD2 mixing ratio (ppb)
#    column 30: gas HCOOH mixing ratio (ppb)
#    column 31: gas RCOOH mixing ratio (ppb)
#    column 32: gas C2O3 mixing ratio (ppb)
#    column 33: gas PAN mixing ratio (ppb)
#    column 34: gas ARO1 mixing ratio (ppb)
#    column 35: gas ARO2 mixing ratio (ppb)
#    column 36: gas ALK1 mixing ratio (ppb)
#    column 37: gas OLE1 mixing ratio (ppb)
#    column 38: gas API1 mixing ratio (ppb)
#    column 39: gas API2 mixing ratio (ppb)
#    column 40: gas LIM1 mixing ratio (ppb)
#    column 41: gas LIM2 mixing ratio (ppb)
#    column 42: gas PAR mixing ratio (ppb)
#    column 43: gas AONE mixing ratio (ppb)
#    column 44: gas MGLY mixing ratio (ppb)
#    column 45: gas ETH mixing ratio (ppb)
#    column 46: gas OLET mixing ratio (ppb)
#    column 47: gas OLEI mixing ratio (ppb)
#    column 48: gas TOL mixing ratio (ppb)
#    column 49: gas XYL mixing ratio (ppb)
#    column 50: gas CRES mixing ratio (ppb)
#    column 51: gas TO2 mixing ratio (ppb)
#    column 52: gas CRO mixing ratio (ppb)
#    column 53: gas OPEN mixing ratio (ppb)
#    column 54: gas ONIT mixing ratio (ppb)
#    column 55: gas ROOH mixing ratio (ppb)
#    column 56: gas RO2 mixing ratio (ppb)
#    column 57: gas ANO2 mixing ratio (ppb)
#    column 58: gas NAP mixing ratio (ppb)
#    column 59: gas XO2 mixing ratio (ppb)
#    column 60: gas XPAR mixing ratio (ppb)
#    column 61: gas ISOP mixing ratio (ppb)
#    column 62: gas ISOPRD mixing ratio (ppb)
#    column 63: gas ISOPP mixing ratio (ppb)
#    column 64: gas ISOPN mixing ratio (ppb)
#    column 65: gas ISOPO2 mixing ratio (ppb)
#    column 66: gas API mixing ratio (ppb)
#    column 67: gas LIM mixing ratio (ppb)
#    column 68: gas DMS mixing ratio (ppb)
#    column 69: gas MSA mixing ratio (ppb)
#    column 70: gas DMSO mixing ratio (ppb)
#    column 71: gas DMSO2 mixing ratio (ppb)
#    column 72: gas CH3SO2H mixing ratio (ppb)
#    column 73: gas CH3SCH2OO mixing ratio (ppb)
#    column 74: gas CH3SO2 mixing ratio (ppb)
#    column 75: gas CH3SO3 mixing ratio (ppb)
#    column 76: gas CH3SO2OO mixing ratio (ppb)
#    column 77: gas CH3SO2CH2OO mixing ratio (ppb)
#    column 78: gas SULFHOX mixing ratio (ppb)

#plot "out/urban_plume_wh_0001_gas.txt" using ($1/3600):12 axes x1y1 with lines title "O3", \
#     "out/urban_plume_wh_0001_gas.txt" using ($1/3600):7 axes x1y1 with lines title "NO2", \
#     "out/urban_plume_wh_0001_gas.txt" using ($1/3600):9 axes x1y1 with lines title "N2O5"

plot "out/urban_plume_nh_0001_gas.txt" using ($1/3600):9 axes x1y1 with lines title "nh N2O5", \
     "out/urban_plume_none_0001_gas.txt" using ($1/3600):9 axes x1y1 with lines title "none N2O5", \
     "out/urban_plume_comp_0001_gas.txt" using ($1/3600):9 axes x1y1 with lines title "comp N2O5"

plot "out/urban_plume_nh_0001_gas.txt" using ($1/3600):7 axes x1y1 with lines title "nh NO2", \
     "out/urban_plume_none_0001_gas.txt" using ($1/3600):7 axes x1y1 with lines title "none NO2", \
     "out/urban_plume_comp_0001_gas.txt" using ($1/3600):7 axes x1y1 with lines title "comp NO2"

#plot "out/urban_plume_wh_0001_gas.txt" using ($1/3600):3 axes x1y1 with lines title "HNO3", \
#     "out/urban_plume_wh_0001_gas.txt" using ($1/3600):24 axes x1y1 with lines title "HCHO", \
#     "out/urban_plume_wh_0001_gas.txt" using ($1/3600):19 axes x1y1 with lines title "SO2", \
#     "out/urban_plume_wh_0001_gas.txt" using ($1/3600):5 axes x1y1 with lines title "NH3"

unset multiplot
