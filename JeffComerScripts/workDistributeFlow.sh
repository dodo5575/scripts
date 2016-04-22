#!/bin/bash

for f in $@; do
    name=${f%.*}

    cmd="browntownf $f $name 0 > $name.log & disown"
    echo "JOB $! $name"
    echo $cmd
    eval $cmd
done
