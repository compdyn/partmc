#!/bin/sh

for f in fig_aging_*.py ; do
    echo
    echo "\vspace{1em}"
    echo "\begin{verbatim}"
    echo $f
    ./$f
    echo "\end{verbatim}"
done
