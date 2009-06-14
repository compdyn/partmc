# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "k"
set ylabel "prob(k)"

set key top right

plot "out/poisson_1_exact.dat" with lines title "mean 1, exact"
replot "out/poisson_1_approx.dat" with points title "mean 1, sampled approx"

replot "out/poisson_4_exact.dat" with lines title "mean 4, exact"
replot "out/poisson_4_approx.dat" with points title "mean 4, sampled approx"

replot "out/poisson_10_exact.dat" with lines title "mean 10, exact"
replot "out/poisson_10_approx.dat" with points title "mean 10, sampled approx"

replot "out/poisson_20_exact.dat" with lines title "mean 20, exact"
replot "out/poisson_20_approx.dat" with points title "mean 20, sampled approx"

replot "out/poisson_30_exact.dat" with lines title "mean 30, exact"
replot "out/poisson_30_approx.dat" with points title "mean 30, sampled approx"
