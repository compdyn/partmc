# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "time / s"
set ylabel "temperature / K"
set y2label "relative humidity"

set key center right

set ytics nomirror
set y2tics

plot "out_amanda_we_eq/condense_flat_0001_env.txt" using 1:2 axes x1y1 with lines title "w.e. temperature", \
     "out_amanda_we_eq/condense_flat_0001_env.txt" using 1:3 axes x1y2 with lines title "w.e. relative humidity", \
     "out_amanda_ne_eq/condense_flat_0001_env.txt" using 1:2 axes x1y1 with lines title "n.e. temperature", \
     "out_amanda_ne_eq/condense_flat_0001_env.txt" using 1:3 axes x1y2 with lines title "n.e. relative humidity"


#plot "out_amanda_eq/up2_hour01_eq_0001_env.txt" using 1:2 axes x1y1 with lines title "w.e. temperature", \
#     "out_amanda_eq/up2_hour01_eq_0001_env.txt" using 1:3 axes x1y2 with lines title "w.e. relative humidity", \

#plot "out_amanda_we/condense_flat_0001_env.txt" using 1:2 axes x1y1 with lines title "w.e. temperature", \
#     "out_amanda_we/condense_flat_0001_env.txt" using 1:3 axes x1y2 with lines title "w.e. relative humidity", \
#     "out_amanda_ne/condense_flat_0001_env.txt" using 1:2 axes x1y1 with lines title "n.e. temperature", \
#     "out_amanda_ne/condense_flat_0001_env.txt" using 1:3 axes x1y2 with lines title "n.e. relative humidity"


#plot "out_amanda_we/condense_0001_env.txt" using 1:2 axes x1y1 with lines title "w.e. temperature", \
#     "out_amanda_we/condense_0001_env.txt" using 1:3 axes x1y2 with lines title "w.e. relative humidity", \
#     "out_amanda_ne/condense_0001_env.txt" using 1:2 axes x1y1 with lines title "n.e. temperature", \
#     "out_amanda_ne/condense_0001_env.txt" using 1:3 axes x1y2 with lines title "n.e. relative humidity"

#plot "out/up2_hour20_eq_0001_env.txt" using 1:2 axes x1y1 with lines title "w.e. temperature", \
#     "out/up2_hour20_eq_0001_env.txt" using 1:3 axes x1y2 with lines title "w.e. relative humidity"

#plot "out_we/out_20/condense_0001_env.txt" using 1:2 axes x1y1 with lines title "w.e. temperature", \
#     "out_we/out_20/condense_0001_env.txt" using 1:3 axes x1y2 with lines title "w.e. relative humidity", \
#     "out_ne/out_20/condense_0001_env.txt" using 1:2 axes x1y1 with lines title "n.e. temperature", \
#     "out_ne/out_20/condense_0001_env.txt" using 1:3 axes x1y2 with lines title "n.e. relative humidity", \
#     "out/up2_hour20_eq_0001_env.txt" using 1:2 axes x1y1 with lines title "w.e. temperature", \
#     "out/up2_hour20_eq_0001_env.txt" using 1:3 axes x1y2 with lines title "w.e. relative humidity"
