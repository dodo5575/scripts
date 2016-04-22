# Author: jcomer2@uiuc.edu

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

proc compute {name moiety structPrefix dcd dcdFreq outDir startFrame} {
    set displayPeriod 200
    set stride 1
    set timestep 1.0
    set selText "name POT"
    set pre "density"

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc
    
    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Read the system size from the xsc file.
    # Note: This only works for lattice vectors along the axes!
    set in [open $xsc r]
    foreach line [split [read $in] "\n"] {
	if {![string match "#*" $line]} {
	    set param [split $line]
	    puts $param
	    set lx [lindex $param 1]
	    set ly [lindex $param 5]
	    set lz [lindex $param 9]
	    break
	}
    }
    puts "NOTE: The system size is $lx $ly $lz.\n"
    close $in
    set origin [list [expr -0.5*$lx] [expr -0.5*$ly] [expr -0.5*$lz]]
    set destin [list [expr 0.5*$lx] [expr 0.5*$ly] [expr 0.5*$lz]]

    # Open the output files.
    set out [open "${outDir}/index_${name}.txt" w]

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcd {
	# Load the trajectory.
	animate delete all
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." $nFrames]

	set outName ${pre}_[trimExtension [trimPath $dcdFile]].dx
	set outFile $outDir/$outName
	volmap density $sel -minmax [list $origin $destin] -checkpoint 0 -o $outFile -allframes -combine avg

	puts $out "$outName $nFrames"

	set nFrames0 [expr $nFrames+$nFrames0]
    }

    $sel delete
    close $out
    mol delete top
}
