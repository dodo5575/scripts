# jcomer2@illinois.edu
package require pmepot

if {$argc < 5} {
    puts "$argv0 name moiety structPrefix outputDir dcdFreq stride dcdFile0 [dcdFile1...]"
    exit
}
set name [lindex $argv 0]
set moiety [lindex $argv 1]
set structPrefix [lindex $argv 2]
set outDir [lindex $argv 3]
set dcdFreq [lindex $argv 4]
set stride [lindex $argv 5]
set dcdList [lrange $argv 6 end]

proc compute {name moiety structPrefix dcd dcdFreq outDir stride} {
    set displayPeriod 200
    set timestep 1.0
    set outFile "$outDir/poten_${name}.dx"
    set pmeSize {128 128 256}

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc pme.xsc

    # Load the system.
    mol load psf $psf pdb $pdb
 
    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcd {
	# Load the trajectory.
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]
    }

    # Run pmePot.
    pmepot -mol top -xscfile $xsc -ewaldfactor 0.25 -grid $pmeSize -dxfile $outFile -frames all
    mol delete top
}

compute $name $moiety $structPrefix $dcdList $dcdFreq $outDir $stride
exit
