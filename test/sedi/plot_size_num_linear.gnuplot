# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "radius / m"
set ylabel "number concentration / (1/m^3)"

set xrange [1e-7:1e-2]
set yrange [0:1.2e9]

plot "out/sedi_part_size_num.txt" using 1:2 title "particle t = 0 hours"
replot "out/sedi_part_size_num.txt" using 1:3 title "particle t = 5 minutes"
replot "out/sedi_part_size_num.txt" using 1:4 title "particle t = 10 minutes"

replot "out/sedi_sect_size_num.txt" using 1:2 with lines title "sectional t = 0 hours"
replot "out/sedi_sect_size_num.txt" using 1:3 with lines title "sectional t = 5 minutes"
replot "out/sedi_sect_size_num.txt" using 1:4 with lines title "sectional t = 10 minutes"
