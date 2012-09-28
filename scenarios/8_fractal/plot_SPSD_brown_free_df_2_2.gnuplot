
# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "dimensionless volume"
set ylabel "dimensionless number concentration"

set xrange [1e-3:10]
set yrange [0:2]

plot "out/part_brown_free_df_2_2_0001_self_preserve.txt" using 1:2 title "Df = 2.2, PartMC", \
     "ref_free_df_2_2_self_preserve_regrid.txt" using 1:2 with lines title "Df = 2.2, Ref"
