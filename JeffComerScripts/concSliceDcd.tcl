# Calculate the current for a trajectory.
# to use: vmd -dispdev text -e concSliceDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set stride 1
set bufK 4
set bufCl 4
set excludeK "(not within $bufK of resname SIN)"
set excludeCl "(not within $bufCl of resname SIN)"
set selTextList [list "name OH2 and $excludeK" "name POT and $excludeK" "name CLA and $excludeCl"]
set selNameList [list "H2O" "K" "Cl"]
set timestep 1.0
set nSegments 30
set binZStart -150.0
set binZEnd 150.0

# Input:
set dir0 /Scr/nanopore/jcomer/myhairpin/nonstick_tail_dcd
set pdb hairpin_E4V_2ns.pdb
set psf pore+dna_E4V.psf
set dcdPrefix $dir0/ecoil_first_4Va
set dcdSuffix ".dcd"
set dcdSet {5 6 7 8 9 10}
# Output:
set outFilePrefix conc_slice_coil_first_4Va[lindex $dcdSet 0]-[lindex $dcdSet end]

set dcdEnd [llength $dcdSet]
# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

# Prepare the segments.
set binDz [expr ($binZEnd-$binZStart)/$nSegments]
set binZ0 {}
set binZ1 {}
for {set i 0} {$i < $nSegments} {incr i} {
	lappend binZ0 [expr $binZStart+$i*$binDz]
	lappend binZ1 [expr $binZStart+($i+1)*$binDz]
}

# Write the bin geometry.
set outBin [open ${outFilePrefix}.bin w]
foreach z0 $binZ0 z1 $binZ1 {
	puts $outBin "$z0 $z1"
}
close $outBin

# Prepare the file names.
set fileNameList {}
foreach name $selNameList {
	lappend fileNameList ${outFilePrefix}.${name}
}

proc computeSegQuantities {selText binZ0 binZ1} {
	set quant {}
	
	# Loop over the segments, computing values for each.
	foreach z0 $binZ0 z1 $binZ1 {
		set sel [atomselect top "z >= $z0 and z < $z1 and (${selText})"]
		lappend quant [$sel num]
		$sel delete
	}
		
	return $quant
}

# Open the output files.
set out {}
foreach name $fileNameList {
	lappend out [open $name w]
}
set outTime [open ${outFilePrefix}.t w]

# Load the system.
mol load psf $psf pdb $pdb

# Loop over the dcd files.
set nFrames0 0
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {
    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Start at frame 0 and move forward, computing
    # the quantities at each step.
    for {set f 0} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
		
	set t [expr ($nFrames0+$f)*$dt]
	puts $outTime $t
	
	# Loop over the selections.
	foreach s $selTextList o $out {
		puts $o [computeSegQuantities $s $binZ0 $binZ1] 
	}
				
	# Show the progress.
	puts "FRAME $f: $t ns"
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

foreach o $out {
	close $o
}
close $outTime
exit



