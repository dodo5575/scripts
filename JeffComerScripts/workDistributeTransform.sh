#!/bin/bash

for f in $@; do
    name=${f##*/}
    out=${name%.*}

    isPot=`echo $out | grep pot`

    if [ -n "$isPot" ]; then
	cou=6.156964
    else
	cou=-6.156964
    fi

    cmd="gridTransformAddCoulomb $f $cou grids/$out.dx"
    echo ""
    echo $cmd
    eval $cmd
done
