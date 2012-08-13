# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / s"
set ylabel "aerosol water mass mixing ratio / (kg/kg)"

set key bottom right

set yrange [0:*]

plot "out/condense_0001_env_aero_time.txt" using 1:($28/($4/(287*$2))) with lines title "aerosol water"
