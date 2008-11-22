#!/bin/sh

for f in fig_*.py ; do
    echo ./$f
    ./$f
done
