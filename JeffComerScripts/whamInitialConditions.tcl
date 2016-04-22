# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set keyList {14 14.5 15 15.5 16 16.5 17 17.5 18 18.5 19}
set selText "resname DMMP"
set occupancy 1.0
set stride 1
set maxShift 0.2
# Input:
set psf DMMP_sol.psf
set dcd output/chemx_dmmp_f110_smd0.dcd
# Output:
set outPrefix wham_dmmp_z

# Load the system.
mol load psf $psf
mol addfile $dcd waitfor all
set nFrames [molinfo top get numframes]
puts [format "Reading %i frames." $nFrames]

# Make the selections.
set selAll [atomselect top all]
set sel [atomselect top $selText]
puts "Performing WHAM on [$sel num] atoms."

# Write the smd files.    
foreach key $keyList {
    # Find the frame with conditions closest to the key state.
    molinfo top set frame 0
    set closeFrame 0
    set closeKey [lindex [measure center $sel weight mass] 2]
    for {set f 1} {$f < $nFrames} {incr f $stride} {
	molinfo top set frame $f
	
	set newKey [lindex [measure center $sel weight mass] 2]
	if {abs($newKey-$key) < abs($closeKey-$key)} {
	    set closeFrame $f
	    set closeKey $newKey
	}
    }
   
    # Choose the closest frame.
    molinfo top set frame $closeKey
    
    $selAll set occupancy 0.0
    $sel set occupancy $occupancy

    set d [expr $key-$closeKey]
    if {abs($d) > $maxShift} {
	puts "Original position: $closeKey"
	puts "Modified position: $key"
	puts "Warning! Moving molecule $d angstroms!"
    }
    $sel moveby [list 0 0 $d]

    # Write the WHAM file.
    $selAll writepdb ${outPrefix}${key}.pdb
}

$sel delete
$selAll delete

mol delete top
exit



