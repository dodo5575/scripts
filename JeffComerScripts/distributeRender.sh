#!/usr/bin/bsh

# Run tachyon jobs, alternating 4 core machines.

nMachines=4
machines=(tbgl-work2 tbgl-work10 tbgl-work8 tbgl-work9)
user=jcomer
workDir=/projects/jcomer/streptavidin/neutravidin/render

mach=0
for file in $@; do
    cmd="cd $workDir; tachyon -add_skylight 1.3 -rescale_lights 0.6 -res 624 11311 -aasamples 3 $file -o $file.tga > $file.log & disown;"
    sshCmd="ssh -l $user ${machines[$mach]} \"$cmd\""
    echo $sshCmd
    eval $sshCmd

    # Increment the machine.
    ((mach++))
    if [ $mach -eq $nMachines ]; then
	mach=0
    fi
done
