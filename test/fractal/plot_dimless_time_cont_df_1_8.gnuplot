
# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set logscale y
set xlabel "dimensionless time"
set ylabel "normalized number concentration"

set xrange [0.01:1000]
set yrange [1e-4:1]

plot "out_dimless_t/part_brown_cont_df_1_8_0001_dimless_t_series.txt" using 1:2 title "Df = 1.8, PartMC", \
     "ref_cont_df_1_8_dimless_time_regrid.txt" using 1:2 with lines title "Df = 1.8, Ref"
