# Calculate the number of charge carriers and current (nA) with multiple
# selections for a trajectory.
# Output format is t(ns) nSel0 ISel0(nA) nSel1 ISel1(nA) ... nTotal ITotal(nA).
# to use: vmd -dispdev text -e multicurrentDCD.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 1000
set radN_K  [expr 3.2]
set radN_Cl [expr 4.2]

set nearSinK "(within $radN_K of resname SIN)"
set nearSinCl "(within $radN_Cl of resname SIN)"
set seg0 "(z >= -127 and z < -102)"
set seg1 "(z >= -102 and z < -25)"
set seg2 "(z >= -25 and z < 0)"
set seg3 "(z >= 0 and z < 25)"
set seg4 "(z >= 25 and z < 102)"
set seg5 "(z >= 102 and z < 127)"
set seg6 "(z < -127 or z >= 127)"

set selSinK0 "name POT and $nearSinK and $seg0"
set selBulkK0 "name POT and not $nearSinK and $seg0"
set selSinK1 "name POT and $nearSinK and $seg1"
set selBulkK1 "name POT and not $nearSinK and $seg1"
set selSinK2 "name POT and $nearSinK and $seg2"
set selBulkK2 "name POT and not $nearSinK and $seg2"
set selSinK3 "name POT and $nearSinK and $seg3"
set selBulkK3 "name POT and not $nearSinK and $seg3"
set selSinK4 "name POT and $nearSinK and $seg4"
set selBulkK4 "name POT and not $nearSinK and $seg4"
set selSinK5 "name POT and $nearSinK and $seg5"
set selBulkK5 "name POT and not $nearSinK and $seg5"
set selSinK6 "name POT and $nearSinK and $seg6"
set selBulkK6 "name POT and not $nearSinK and $seg6"

set selSinCl0 "name CLA and $nearSinCl and $seg0"
set selBulkCl0 "name CLA and not $nearSinCl and $seg0"
set selSinCl1 "name CLA and $nearSinCl and $seg1"
set selBulkCl1 "name CLA and not $nearSinCl and $seg1"
set selSinCl2 "name CLA and $nearSinCl and $seg2"
set selBulkCl2 "name CLA and not $nearSinCl and $seg2"
set selSinCl3 "name CLA and $nearSinCl and $seg3"
set selBulkCl3 "name CLA and not $nearSinCl and $seg3"
set selSinCl4 "name CLA and $nearSinCl and $seg4"
set selBulkCl4 "name CLA and not $nearSinCl and $seg4"
set selSinCl5 "name CLA and $nearSinCl and $seg5"
set selBulkCl5 "name CLA and not $nearSinCl and $seg5"
set selSinCl6 "name CLA and $nearSinCl and $seg6"
set selBulkCl6 "name CLA and not $nearSinCl and $seg6"

set selText [list $selSinK0 $selSinK1 $selSinK2 $selSinK3 $selSinK4 $selSinK5 $selSinK6 $selBulkK0 $selBulkK1 $selBulkK2 $selBulkK3 $selBulkK4 $selBulkK5 $selBulkK6]
lappend selText $selSinCl0 $selSinCl1 $selSinCl2 $selSinCl3 $selSinCl4 $selSinCl5 $selSinCl6 $selBulkCl0 $selBulkCl1 $selBulkCl2 $selBulkCl3 $selBulkCl4 $selBulkCl5 $selBulkCl6

set charge [list 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0 -1.0]
set timestep 1.0

# Input:
set pdb  poreEmpty-all.pdb
set psf  poreEmpty-all.psf
set dcdPrefix  /Scr/alek2/20_DS4/E4V-DS4-noDNA-
set dcdSuffix ".dcd"
set dcdSet {1 2 3 4}
set xsc E4V-DS4-noDNA-4.restart.xsc
# Output:
set outFile seg_hairpin_empty.txt

set dcdEnd [llength $dcdSet]
# Get the time change between frames in femtoseconds.
set dt [expr $timestep*$dcdFreq]

# Return the current and number of charge carriers in the selection for the frame.
proc computeCurrent {frameCurr selText charge lz dt} {
	set frameLast [expr $frameCurr-1]
	
	# Get the current position of the selection.
	molinfo top set frame $frameCurr
	set sel [atomselect top $selText]
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
		set currentZ [expr $currentZ + $charge*$dz/($lz*$dt)*1.60217733e5]
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

# Open the output file.
set out [open $outFile w]
#puts $out "sum of q*v for $psf with trajectory $dcd"
#puts $out "t(ns) I(A)"

# Load the system.
mol load psf $psf pdb $pdb

# Loop over the dcd files.
set nFrames0 0
for {set dcd 0} {$dcd < $dcdEnd} {incr dcd} {

# Load the trajectory.
animate delete all
set dcdNum [lindex $dcdSet $dcd]
mol addfile "${dcdPrefix}${dcdNum}${dcdSuffix}" waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Move forward, computing
# current at each step.
for {set f 1} {$f < $nFrames} {incr f} {
	set fLast [expr $f-1]
	
	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f+0.5)*$dt*1.e-6]
	
	# Write the time.
	puts -nonewline $out $t
	puts -nonewline "FRAME $f: $t"
		
	# Compute the current for each selection.		
	set currentZTotal 0.0
	set numTotal 0
	foreach s $selText q $charge {
		#Write the number of carriers and the current for this selection.
		set data [computeCurrent $f $s $q $lz $dt]
		set num [lindex $data 0]
		set currentZ [lindex $data 1]
		puts -nonewline $out " $num $currentZ"
		puts -nonewline " $currentZ"
		
		set numTotal [expr $numTotal + $num]
		set currentZTotal [expr $currentZTotal + $currentZ]
	}
	
	# Write the total number of carriers and total current.
	puts -nonewline $out " $numTotal"
	puts $out " $currentZTotal"
	puts " $currentZTotal"	
}
set nFrames0 [expr $nFrames+$nFrames0]
}

close $out
exit




