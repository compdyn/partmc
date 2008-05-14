#!/bin/sh

for f in fig_*.py ; do
    echo ./$f $1
    ./$f $1
done
