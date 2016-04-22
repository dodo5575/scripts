# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set keyList {14.5 15 15.5 16 16.5 17 17.5 18 18.5 19}
set latLength 25
set latN 5
set selText "resname DMMP"
set occupancy 1.0
set stride 1
set maxShift 0.2
# Input:
set psf DMMP_sol.psf
set pdb DMMP_sol.pdb
set dcdDir output
set dcdPrefix chemx_dmmp_f
set dcdSet {20 50 60 80 100 110}
set dcdSuffix _smd0.dcd
# Output:
set outPrefix map/dmmp_pos

# Load the system.
mol load psf $psf pdb $pdb

# Determine the node positions.
set nodePos {}
foreach z $keyList {
    for {set ix 0} {$ix < $latN} {incr ix} {
	    for {set iy 0} {$iy < $latN} {incr iy} {
		set x [expr $latLength*(1.0*$ix/$latN-0.5)]
		set y [expr $latLength*(1.0*$iy/$latN-0.5)]
		lappend nodePos [list $x $y $z]
	    }
    }
}
#puts $nodePos

# Make the selections.
set selAll [atomselect top all]
set sel [atomselect top $selText]
puts "Performing WHAM on [$sel num] atoms."

#set closeFile ${dcdPrefix}[lindex $dcdSet 0]${dcdSuffix}
# Set the initial search conditions.
set n [llength $nodePos]
for {set i 0} {$i < $n} {incr i} {
    set closeFrame($i) 0
    set closeDist($i) 1e20
}

# Load all the frames.
animate delete all
foreach dcdNum $dcdSet {
    # Load the trajectory.
    set dcd $dcdDir/${dcdPrefix}${dcdNum}${dcdSuffix}
    mol addfile $dcd waitfor all
}
set nFrames [molinfo top get numframes]
puts [format "Read %i frames." $nFrames]

# Loop through the frames, looking for the closest matches.
for {set f 0} {$f < $nFrames} {incr f $stride} {
    molinfo top set frame $f
    
    set currPos [measure center $sel weight mass]
    for {set i 0} {$i < $n} {incr i} {
	set pos [lindex $nodePos $i]
	set d [veclength [vecsub $currPos $pos]]
	
	if {$d < $closeDist($i)} {
	    set closeDist($i) $d
	    set closeFrame($i) $f
	}
    }
}

# Write the smd files.    
set out [open $outPrefix.txt w]
for {set i 0} {$i < $n} {incr i} {
    set pos [lindex $nodePos $i]
    
    molinfo top set frame $closeFrame($i)
    if {$closeDist($i) > $maxShift} {
	puts "Warning! Moving molecule $closeDist($i) angstroms for position $pos!"
    }
    puts "System $i, position $pos: moving $closeDist($i) angstroms."
    
    $selAll set occupancy 0.0
    $sel set occupancy $occupancy
    set currPos [measure center $sel weight mass]
    $sel moveby [vecsub $pos $currPos]

    # Write the WHAM file.
    $selAll writepdb ${outPrefix}${i}.pdb
    
    puts $out "$i $pos $closeDist($i)"
}
close $out

$sel delete
$selAll delete
mol delete top
exit



