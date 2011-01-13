# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "diameter / m"
set ylabel "number concentration / (1/m^3)"

set xrange [1e-9:1e-6]

set key top left

plot "out/parallel_mix_0001_0001_aero_size_num.txt" using 1:26 with lines title "proc 1", \
     "out/parallel_mix_0001_0002_aero_size_num.txt" using 1:26 with lines title "proc 2", \
     "out/parallel_mix_0001_0003_aero_size_num.txt" using 1:26 with lines title "proc 3", \
     "out/parallel_mix_0001_0004_aero_size_num.txt" using 1:26 with lines title "proc 4", \
     "out/sect_aero_size_num.txt" using 1:26 with lines linewidth 5 title "sect"
