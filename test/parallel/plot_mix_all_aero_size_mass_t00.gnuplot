# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale
set xlabel "diameter / m"
set ylabel "mass concentration / (kg/m^3)"

set xrange [1e-9:1e-6]
set yrange [1e-14:1e-6]

set key top left

plot "out/parallel_mix_0001_0001_aero_size_mass.txt" using 1:2 with lines title "proc 01", \
     "out/parallel_mix_0001_0002_aero_size_mass.txt" using 1:2 with lines title "proc 02", \
     "out/parallel_mix_0001_0003_aero_size_mass.txt" using 1:2 with lines title "proc 03", \
     "out/parallel_mix_0001_0004_aero_size_mass.txt" using 1:2 with lines title "proc 04", \
     "out/parallel_mix_0001_0005_aero_size_mass.txt" using 1:2 with lines title "proc 05", \
     "out/parallel_mix_0001_0006_aero_size_mass.txt" using 1:2 with lines title "proc 06", \
     "out/parallel_mix_0001_0007_aero_size_mass.txt" using 1:2 with lines title "proc 07", \
     "out/parallel_mix_0001_0008_aero_size_mass.txt" using 1:2 with lines title "proc 08", \
     "out/parallel_mix_0001_0009_aero_size_mass.txt" using 1:2 with lines title "proc 09", \
     "out/parallel_mix_0001_0010_aero_size_mass.txt" using 1:2 with lines title "proc 10", \
     "out/sect_aero_size_mass.txt" using 1:2 with lines linewidth 5 title "sect"
