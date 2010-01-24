# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time (s)"
set ylabel "aerosol number concentration (particles/m^3)"
set key right bottom

plot "out/aero_total.txt" using 1:2 title "particle", \
     "out/nucleate_ode_aero_number.txt" using 1:2 w l title "ODE"
