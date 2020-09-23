# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale x
set xlabel "diameter / m"
set ylabel "number concentration / (1/m^3)"

set xrange [1e-7:1e-1]
set yrange [0:1.2e9]

plot "out/sedi_part_0001_aero_size_num.txt" using 1:2 title "particle t = 0 hours", \
     "out/sedi_part_0001_aero_size_num.txt" using 1:3 title "particle t = 5 minutes", \
     "out/sedi_part_0001_aero_size_num.txt" using 1:4 title "particle t = 10 minutes", \
     "out/sedi_sect_aero_size_num.txt" using 1:2 with lines title "sectional t = 0 hours", \
     "out/sedi_sect_aero_size_num.txt" using 1:3 with lines title "sectional t = 5 minutes", \
     "out/sedi_sect_aero_size_num.txt" using 1:4 with lines title "sectional t = 10 minutes"
