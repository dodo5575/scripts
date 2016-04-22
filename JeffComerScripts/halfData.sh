#!/bin/bash

for file in $@; do
	awk 'NR % 2 == 0 {print $0}' $file > tmp0.dat
	awk 'NR % 2 == 0 {print $0}' $file > tmp1.dat
	paste tmp0.dat tmp1.dat | awk '{print 0.5*($1+$4),($2+$5),0.5*($3+$6)}' > $file.comb
	echo $file.comb
done
