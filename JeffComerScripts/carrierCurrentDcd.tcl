# Calculate the current (nA) with multiple
# selections for a trajectory (ns).
# to use: vmd -dispdev text -e carrierCurrentDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set selText [list "name POT" "name CLA"]
set charge [list 1.0 -1.0]
set name [list "K" "Cl"]
set timestep 1.0
set stride 1

# Input:
set psf open_1.5V.psf
set pdb open_1.5V.pdb
set dcdPrefix /Scr/nanopore/jcomer/myhairpin/open_1.5V_dcd/open_1.5V
set dcdSuffix .dcd
set dcdSet {0 1 2}
set xsc output/eq4.restart.xsc
# Output:
set outFileSuffix open_1.5V_[lindex $dcdSet 0]-[lindex $dcdSet end].dat

set dcdEnd [llength $dcdSet]
# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

# Return the current and number of charge carriers in the selection for the frame.
proc computeCurrent {frameCurr sel charge lz dt} {
	set frameLast [expr $frameCurr-1]
	
	# Get the current position of the selection.
	molinfo top set frame $frameCurr
	set z1 [$sel get z]
		
	# Get the last position of the selection.
	molinfo top set frame $frameLast
	set z0 [$sel get z]
		
	# Get the number of charge carriers.
	set num [$sel num]
	
	# Find the displacements in the z-direction and compute the current.
	set currentZ 0.0
	foreach a0 $z0 a1 $z1 {
		# Compensate for jumps across the periodic cell.
		set dz [expr $a1-$a0]
		if {[expr $dz > 0.5*$lz]} {set dz [expr $dz-$lz]}
		if {[expr $dz <-0.5*$lz]} {set dz [expr $dz+$lz]}

		# Compute the current in nanoamperes.
		set currentZ [expr $currentZ + $charge*$dz/($lz*$dt)*1.60217733e-1]
	}
	
	return [list $num $currentZ]
}

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

# Open the output files.
set out {}
foreach n $name {
	lappend out [open "curr${n}_${outFileSuffix}" w]
}
set outTotal [open "curr_${outFileSuffix}" w]

# Load the system.
mol load psf $psf pdb $pdb
set sel {}
foreach st $selText {
	lappend sel [atomselect top $st]
}

# Loop over the dcd files.
set nFrames0 0
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

    # Load the trajectory.
    animate delete all
    set dcdNum [lindex $dcdSet $dcd]
    mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Move forward, computing
    # current at each step.
    for {set f 1} {$f < $nFrames} {incr f} {
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f+0.5)*$dt]
	
	# Write the time.
	puts -nonewline "FRAME $f: $t"
		
	# Compute the current for each selection.		
	set currentZTotal 0.0
	foreach s $sel q $charge o $out {
	    #Write the number of carriers and the current for this selection.
	    set data [computeCurrent $f $s $q $lz $dt]
	    set currentZ [lindex $data 1]
	    puts $o "$t $currentZ"
	    puts -nonewline " $currentZ"
		
	    set currentZTotal [expr $currentZTotal + $currentZ]
	}
	
	# Write the total current.
	puts $outTotal "$t $currentZTotal"
	puts " $currentZTotal"	
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

foreach s $sel o $out {
	$s delete
	close $o
}
exit




