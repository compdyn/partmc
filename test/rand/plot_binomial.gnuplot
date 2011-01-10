# run from inside gnuplot with:
# load "<filename>.gnuplot"
# or from the commandline with:
# gnuplot -persist <filename>.gnuplot

set xlabel "k"
set ylabel "prob(k)"

set key top right

plot "out/binomial_1_exact.dat" with lines title "n = 10, p = 0.3, exact", \
     "out/binomial_1_approx.dat" with points title "n = 10, p = 0.3, sampled approx", \
     "out/binomial_2_exact.dat" with lines title "n = 10, p = 0.7, exact", \
     "out/binomial_2_approx.dat" with points title "n = 10, p = 0.7, sampled approx", \
     "out/binomial_3_exact.dat" with lines title "n = 30, p = 0.5, exact", \
     "out/binomial_3_approx.dat" with points title "n = 30, p = 0.5, sampled approx"
