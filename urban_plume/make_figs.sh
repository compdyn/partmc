#!/bin/sh

for f in fig_*.py ; do
    echo
    echo "\vspace{1em}"
    echo "\begin{verbatim}"
    echo $f
    ./$f
    echo "\end{verbatim}"
done
