# Calculate the number of charge carriers and current (nA) with multiple
# selections for a trajectory.
# Output format is t(ns) nSel0 ISel0(nA) nSel1 ISel1(nA) ... nTotal ITotal(nA).
# to use: vmd -dispdev text -e dnaMulticurrentDCD.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set dcdFreq 5000
set radN_K  3.76125
set radN_Cl 4.2675

set nearDnaK "(within $radN_K of segname ADNA BDNA)"
set nearSinK "(within $radN_K of resname SIN)"
set nearDnaCl "(within $radN_Cl of segname ADNA BDNA)"
set nearSinCl "(within $radN_Cl of resname SIN)"

set selBothK "name POT and $nearSinK and $nearDnaK"
set selDnaK "name POT and not $nearSinK and $nearDnaK"
set selSinK "name POT and $nearSinK and not $nearDnaK"
set selBulkK "name POT and not $nearSinK and not $nearDnaK"

set selBothCl "name CLA and $nearSinCl and $nearDnaCl"
set selDnaCl "name CLA and not $nearSinCl and $nearDnaCl"
set selSinCl "name CLA and $nearSinCl and not $nearDnaCl"
set selBulkCl "name CLA and not $nearSinCl and not $nearDnaCl"

set selText [list $selBothK $selDnaK $selSinK $selBulkK $selBothCl $selDnaCl $selSinCl $selBulkCl]
set charge [list 1.0 1.0 1.0 1.0 -1.0 -1.0 -1.0 -1.0]
set timestep 1.0

# Input:
set pdb  DS3_pore+dna-all.pdb
set psf  DS3_pore+dna-all.psf
set dcdPrefix  /Scr/alek/DS3-MORE/output/electricDS3-E0.1-
set dcdSuffix ".dcd"
set dcdSet {12 13 14 15}
set xsc /Scr/alek/DS3-MORE/output/electricDS3-E0.1-15.restart.xsc
# Output:
set outFile dna_curr_DS3_E0.1.txt

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




