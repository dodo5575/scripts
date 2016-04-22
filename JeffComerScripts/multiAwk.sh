#!/usr/bin/bsh

name="pot-chl"
last=38
prefix="umbrella_${name}"
suffix=".dat"
outPrefix="three_${name}"

for ((i=0;i<=$last;i+=1)); do
	awk '{print $1,0,0,$2}' ${prefix}${i}${suffix} > ${outPrefix}${i}${suffix}
done
