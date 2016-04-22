# jcomer2@uiuc.edu
proc compute {name moiety structPrefix dcdList dcdFreq outDir stride} {
    set timestep 1.0
    set selText "protein and name CA"
    set outFile fit_${name}.dcd

    # Input:
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb
    set xsc $structPrefix.xsc

    # Load the system.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    set all [atomselect top all]
    puts "Loaded the structure `$psf' `$pdb'."
    animate delete all
 
    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcdList {
	# Load the trajectory.
	mol addfile $dcdFile type dcd step $stride waitfor all
	set nFrames [molinfo top get numframes]
	puts [format "Reading %i frames." [expr {$nFrames-$nFrames0}]]
	incr nFrames0 $nFrames
    }
    set totalFrames [molinfo top get numframes]
    puts "Loaded [llength $dcdList] dcd files."
    
    # Regenerate each frame.
    for {set f 0} {$f < $totalFrames} {incr f} {
	molinfo top set frame $f
	set cen [measure center $sel weight mass]
	$all moveby [vecinvert $cen]
    }

    animate write dcd $outDir/$outFile beg 0 end -1 sel $all top
    puts "Wrote the output dcd."

    $sel delete
    $all delete
    mol delete top
}
