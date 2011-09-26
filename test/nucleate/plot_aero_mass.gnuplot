# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / s"
set ylabel "aerosol SO4 mass concentration / (kg/m^3)"
set key right bottom

plot "out/nucleate_part_0001_aero_time.txt" using 1:4 title "particle", \
     "out/nucleate_ode_aero_time.txt" using 1:3 w l title "ODE"
