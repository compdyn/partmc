# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set logscale
set xlabel "diameter / m"
set ylabel "mass concentration / (kg/m^3)"

set key left top

set xrange [1e-7:1e-1]
set yrange [1e-13:1e0]

plot "out/sedi_part_0001_aero_size_mass.txt" using 1:2 title "particle t = 0 hours", \
     "out/sedi_part_0001_aero_size_mass.txt" using 1:3 title "particle t = 5 minutes", \
     "out/sedi_part_0001_aero_size_mass.txt" using 1:4 title "particle t = 10 minutes", \
     "out/sedi_sect_aero_size_mass.txt" using 1:2 with lines title "sectional t = 0 hours", \
     "out/sedi_sect_aero_size_mass.txt" using 1:3 with lines title "sectional t = 5 minutes", \
     "out/sedi_sect_aero_size_mass.txt" using 1:4 with lines title "sectional t = 10 minutes"
