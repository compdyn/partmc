# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "diameter / m"
set ylabel "mass concentration / (kg/m^3)"

set xrange [1e-9:1e-6]
#set yrange [0:5e8]

plot "out/average_0001_aero_size_mass.txt" using 1:2 with lines title "base", \
     "out/average_comp_0001_aero_size_mass.txt" using 1:2 with points pointsize 3 title "comp average", \
     "out/average_sizenum_0001_aero_size_mass.txt" using 1:2 with lines title "size(num) average", \
     "out/average_sizevol_0001_aero_size_mass.txt" using 1:2 with lines title "size(vol) average", \
     "out/average_compsizenum_0001_aero_size_mass.txt" using 1:2 with lines title "comp then size(num) average", \
     "out/average_compsizevol_0001_aero_size_mass.txt" using 1:2 with points title "comp then size(vol) average"
