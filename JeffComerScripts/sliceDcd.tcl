# Calculate the current for a trajectory.
# to use: vmd -dispdev text -e sliceDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selTextList [list "name OH2 and resname TIP3" "name POT" "name CLA" "segname HAIR"]
set selNameList [list "H2O" "K" "Cl" "DNA"]
set timestep 1.0
set nSegments 60
set binZStart -150.0
set binZEnd 150.0

# Input:
set pdb  pore+dna-all.pdb
set psf  pore+dna-all.psf
set dcdPrefix  run
set dcdSuffix "_4V.dcd"
set dcdSet {10 11 12}
# Output:
set outFilePrefix slice_conform

set dcdEnd [llength $dcdSet]
# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

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
		set sel [atomselect top "(${selText}) and z >= $z0 and z < $z1"]
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
set nFrames0 1
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

# Load the trajectory.
animate delete all
set dcdNum [lindex $dcdSet $dcd]
mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Start at frame 0 and move forward, computing
# the quantities at each step.
for {set f 0} {$f < $nFrames} {incr f} {
	molinfo top set frame $f
		
	set t [expr ($nFrames0+$f)*$dt*1.e-6]
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




