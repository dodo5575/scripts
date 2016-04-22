# jcomer2@illinois.edu
if {$argc < 5} {
    puts "$argv0 name moiety structPrefix outputDir stride dcdFile0 \[dcdFile1...\]"
    exit
}

set name [lindex $argv 0]
set moiety [lindex $argv 1]
set structPrefix [lindex $argv 2]
set outDir [lindex $argv 3]
set stride [lindex $argv 4]
set dcdList [lrange $argv 5 end]

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

proc compute {name moiety structPrefix dcd outDir stride} {
    #set gridFile best_pmf_no.dx
    set displayPeriod 200
    set selText "not (water or resname SIN SIO2)"
    if {$stride < 1} {set stride 1}

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

compute $name $moiety $structPrefix $dcdList $outDir $stride
exit
