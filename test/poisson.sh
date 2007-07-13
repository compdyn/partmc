#!/bin/sh

cat <<ENDINFO

Poisson Sample Test
-------------------

This test generates samples from the Poisson random number generator
used in the main code. The resulting sampled distribution is plotted
against the analytical expression for the Poisson PDF. The Poisson
generator used at the moment is only approximate, so the two plots
should be similar but not identical. The PDFs are plotted for several
different values of the rate parameter.

ENDINFO
sleep 1

echo ./poisson_sample 1 50 10000000
./poisson_sample 1 50 10000000 > out/poisson_1.d

echo ./poisson_sample 4 50 10000000
./poisson_sample 4 50 10000000 > out/poisson_4.d

echo ./poisson_sample 10 50 10000000
./poisson_sample 10 50 10000000 > out/poisson_10.d

echo ./poisson_sample 20 50 10000000
./poisson_sample 20 50 10000000 > out/poisson_20.d

echo ./poisson_sample 30 50 10000000
./poisson_sample 30 50 10000000 > out/poisson_30.d

echo Plotting Poisson samples
gnuplot -persist <<ENDPOISSON
set xlabel "n"
set ylabel "Prob(n)"
set title "Sampled and exact Poisson PDFs"
plot "out/poisson_1.d" using 1:3 title "rate 1 (sampled)"
replot "out/poisson_1.d" using 1:2 with lines title "rate 1 (exact)"
replot "out/poisson_4.d" using 1:3 title "rate 4 (sampled)"
replot "out/poisson_4.d" using 1:2 with lines title "rate 4 (exact)"
replot "out/poisson_10.d" using 1:3 title "rate 10 (sampled)"
replot "out/poisson_10.d" using 1:2 with lines title "rate 10 (exact)"
replot "out/poisson_20.d" using 1:3 title "rate 20 (sampled)"
replot "out/poisson_20.d" using 1:2 with lines title "rate 20 (exact)"
replot "out/poisson_30.d" using 1:3 title "rate 30 (sampled)"
replot "out/poisson_30.d" using 1:2 with lines title "rate 30 (exact)"
set terminal postscript eps
set output "out/poisson_plot.eps"
replot
ENDPOISSON
epstopdf out/poisson_plot.eps
