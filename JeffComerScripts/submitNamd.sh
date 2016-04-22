#!/usr/bin/bsh

for x in 0 1 2 3; do
	num=$(($1*4 + $x))
	cmd="namd2cvs corr_basepair_null_pos${num}.namd > corr_basepair_null_pos${num}.log & disown"
	echo $cmd
	eval $cmd
done 
