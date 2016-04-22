#!/bin/bash

for f in *.C; do
    prefix=${f%.*}
    [ $f -ot $prefix ] && continue;
    cmd="g++ -O3 -Wall -lgsl -lgslcblas -fopenmp $f -o $prefix"
    echo $cmd
    eval $cmd
done
