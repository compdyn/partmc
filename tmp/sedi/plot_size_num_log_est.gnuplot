# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale
set xlabel "diameter / m"
set ylabel "number concentration / (1/m^3)"

set xrange [1e-7:1e-2]
set yrange [1e-3:1e10]

plot "out/sedi_part_size_num.txt" using 1:2 title "particle t = 0 hours", \
     "out/sedi_part_size_num.txt" using 1:3 title "particle t = 5 minutes", \
     "out/sedi_part_size_num.txt" using 1:4 title "particle t = 10 minutes", \
     "out/sedi_sect_size_num.txt" using 1:2 with lines title "sectional t = 0 hours", \
     "out/sedi_sect_size_num.txt" using 1:3 with lines title "sectional t = 5 minutes", \
     "out/sedi_sect_size_num.txt" using 1:4 with lines title "sectional t = 10 minutes", \
     "out/1_num.txt" using 1:3 with linespoints title "est, t = 0", \
     "out/2_num.txt" using 1:3 with linespoints title "est, t = 5", \
     "out/3_num.txt" using 1:3 with linespoints title "est, t = 10"
