#!/bin/sh

cat README
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

echo ./plot.py
./plot.py
echo Now view out/poisson.pdf
