# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / s"
set ylabel "aerosol number concentration / (1/m^3)"
set key right bottom

plot "out/nucleate_part_0001_aero_time.txt" using 1:2 title "particle", \
     "out/nucleate_ode_aero_time.txt" using 1:2 w l title "ODE"
