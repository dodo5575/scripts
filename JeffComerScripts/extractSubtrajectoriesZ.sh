#!/bin/bash
# Author: Jeff Comer <jcomer2@illinois.edu>

# Output directory
outDir=data
# Length of subtrajectories
subSteps=600

for f in $@; do
    fileName=${f##*/}
    name=${fileName%.*}

    awkCmd="{if (NR % $subSteps == 0) printf \"%.12g\n\", \$5; else printf \"%.12g \", \$5; }"

    #echo $awkCmd
    echo $f
    grep 'SMD[[:space:]][[:space:]]' $f | awk "$awkCmd" > data/$name.z.dat
done
