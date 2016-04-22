#!/bin/bash

for f in *.C; do
    prefix=${f%.*}
    [ ! -f $prefix ] && continue; 
    cmd="rm -f $prefix"
    echo $cmd
    eval $cmd
done
