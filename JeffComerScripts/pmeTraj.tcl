# jcomer2@uiuc.edu
package require pmepot

proc compute {name moiety structPrefix dcd dcdFreq outDir stride} {
    set displayPeriod 200
    set timestep 1.0
    set outFile "$outDir/pot_${name}.dx"
    set pmeSize {96 96 192}

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc

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
