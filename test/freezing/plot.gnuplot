array files = ['01', '02', '03', '04', '05', '06', '07']
a=system("ls -1 *.txt | cut -d_ -f4")
array SUM[|files|]
do for [i=1:|files|] {
  stats 'freezing_part_0001_000000'.files[i].'_aero_particles.txt' using 6 nooutput
  SUM[i] = STATS_sum
}
set style fill solid
set boxwidth 0.5
set xlabel 'output step'
set ylabel '# frozen'
set yrange [0:]
plot SUM using 1:2:xticlabels(files[column(0)+1]) with boxes