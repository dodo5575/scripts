#!/usr/bin/bsh

machines=(tbgl-work2 tbgl-work8 tbgl-work9 tbgl-work10)
name=corr1_basepair_null_pos
user=jcomer
workDir=/projects/jcomer/basepair/correlation_force/namd

for mach in 0 1 2 3; do
	namdCmd="cd $workDir;";
	for x in 0 1 2 3; do
		num=$(($mach*4 + $x))
		cmd="namd2cvs ${name}${num}.namd > ${name}${num}.log & disown;"
		echo $cmd
		namdCmd="$namdCmd $cmd"
	done
	
	sshCmd="ssh -l $user ${machines[$mach]} \"$namdCmd\""
	echo $sshCmd
	eval $sshCmd
done
