# jcomer2@illinois.edu

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

proc compute {name moiety structPrefix dcd dcdFreq outDir startFrame} {
    #set gridFile best_pmf_no.dx
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set selText "not (water or resname SIN SIO2)"
    if {$startFrame < 1} {set startFrame 0}

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
   
     # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    puts "NOTE: Retaining [$sel num] atoms."

    $sel writepsf $outDir/nw_[trimPath $psf]
    $sel writepdb $outDir/nw_[trimPath $pdb]
    
    # Loop over the dcd files.
    foreach dcdFile $dcd {
	animate delete all

	# Load the trajectory.
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	# Write the dcd file.
	set dcdName [trimPath $dcdFile]
	animate write dcd "$outDir/nw_${dcdName}" beg 0 end -1 waitfor all sel $sel top
    }
}
