# Calculate a generalized currrent (1/ns).
# to use: vmd -dispdev text -e generalizedCurrentDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname ADNA and abs(z) < 1.5"
set systemLength 3.0
set chargeQuant mass
set currentFactor 1.0
set dcdFreq 1000
set timestep 1.0
set stride 10

# Input:
set dir /Scr/7day-2/alek/OSCILLATIONS
set psf $dir/polyG25_1M_pore_a30.psf
set pdb $dir/polyG25_1M_pore_a30.pdb
set dcdPrefix $dir/eq
set dcdSuffix a.15V-polyG25_1M.signal_H.dcd
set dcdSet {03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19}
set xsc $dir/eq19a.15V-polyG25_1M.signal_H.xsc
# Output:
set outName polyG25_signal[lindex $dcdSet 0]-[lindex $dcdSet end]

set dcdEnd [llength $dcdSet]
# Get the time change between frames in nanoseconds.
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

# Open the output files.
set out [open mcurr_${outName}.dat w]
set outNum [open num_${outName}.dat w]

# Return the current and number of charge carriers in the selection for the frame.
proc computeCurrent {frameCurr selText chargeQuant lz dt} {
    set frameLast [expr $frameCurr-1]
    
    set sel [atomselect top $selText]
    set charge [$sel get $chargeQuant]

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
    foreach a0 $z0 a1 $z1 q $charge {
	# Compensate for jumps across the periodic cell.
	set dz [expr $a1-$a0]
	if {[expr $dz > 0.5*$lz]} {set dz [expr $dz-$lz]}
	if {[expr $dz <-0.5*$lz]} {set dz [expr $dz+$lz]}

	# Compute the current in nanoamperes.
	set currentZ [expr $currentZ + $q*$dz/($lz*$dt)]
    }

    $sel delete
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

    # Move forward, computing
    # current at each step.
    for {set f 1} {$f < $nFrames} {incr f} {
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f+0.5)*$dt]
	
	# Compute the current for each selection.		
	# Write the number of carriers and the current for this selection.
	set data [computeCurrent $f $selText $chargeQuant $systemLength $dt]
	set num [lindex $data 0]
	set currentZ [expr $currentFactor*[lindex $data 1]]
	puts $out "$t $currentZ"
	puts $outNum "$t $num"
	
	# Write the results.
	puts "FRAME $f: $t $num $currentZ"
		
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
close $outNum
exit




