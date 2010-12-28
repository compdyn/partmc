# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / s"
set ylabel "gas H2SO4 mixing ratio / ppb"

plot "out/gas.txt" using 1:2 title "particle", \
     "out/nucleate_ode_gas.txt" using 1:2 w l title "ODE"
