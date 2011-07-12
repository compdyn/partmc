# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "diameter / m"
set ylabel "mass concentration / (kg/m^3)"

set key left top

set xrange [1e-7:1e-2]
set yrange [0:0.04]

#plot "out/sedi_part_size_mass.txt" using 1:2 title "particle t = 0 hours"
#replot "out/sedi_part_size_mass.txt" using 1:3 title "particle t = 5 minutes"
plot "out/sedi_part_size_mass.txt" using 1:4 title "particle t = 10 minutes"

#replot "out/sedi_sect_size_mass.txt" using 1:2 with lines title "sectional t = 0 hours"
#replot "out/sedi_sect_size_mass.txt" using 1:3 with lines title "sectional t = 5 minutes"
replot "out/sedi_sect_size_mass.txt" using 1:4 with lines title "sectional t = 10 minutes"

#replot "out/1_mass.txt" using 1:3 with linespoints title "est, t = 0"
#replot "out/2_mass.txt" using 1:3 with linespoints title "est, t = 5"
replot "out/3_mass.txt" using 1:3 with linespoints title "est, t = 10"

#replot "out/1_mass.txt" using 1:4 with linespoints title "summed, t = 0"
#replot "out/2_mass.txt" using 1:4 with linespoints title "summed, t = 5"
replot "out/3_mass.txt" using 1:4 with linespoints title "summed, t = 10"
