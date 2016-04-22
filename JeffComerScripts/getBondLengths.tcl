# Extract the coordinates of the chosen atoms.
# Use with: vmd -dispdev text -e getBonds.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set cellA {75.95 0.0 0.0}
set cellB {37.975 65.7746 0.0}
set cellC {0.0 0.0 286.0}
#Input:
set psf pore1.2.psf
set pdb pore1.2.pdb
#Output:
set outFile pore1.2_bond_length.dat

# Determine the basis vectors.
source vector.tcl
set basis [matTranspose [list $cellA $cellB $cellC]]
set basisInv [matInvert $basis]

# Read the bonds from a psf file and return a
# list of zero-based indices.
proc readBonds {psf} {
    set in [open $psf r]
    
    set readingBonds 0
    set bond {}
    foreach line [split [read $in] \n] {
	set tok [concat $line]
	if {[llength $tok] <= 1} {
	    # Quit reading the bonds.
	    set readingBonds 0
	    continue
	}

	if {[string match "*!NBOND*" $line]} {
	    # Start reading bonds.
	    set readingBonds 1
	} elseif {[string match "*!NTHETA*" $line]} {
	    # Quit reading the bonds.
	    set readingBonds 0
	} elseif {$readingBonds} {
	    # Extract a set of bond definitions.
	    foreach {a0 a1} $tok {
		lappend bond [list [expr $a0-1] [expr $a1-1]]
	    }
	}
    }

    close $in
    return $bond
}


# Read the bonds.
foreach zero {0} {set bond [readBonds $psf]}

# Load the molecule.
mol load psf $psf pdb $pdb
set all [atomselect top all]
foreach zero {0} {set pos [$all get {x y z}]}
$all delete

# Map into lattice space.
set latPos {}
foreach r $pos {
    lappend latPos [vecTransform $basisInv $r]
}

# Get the bond lengths in lattice space.
set bondLength {}
foreach b $bond {
    foreach {a0 a1} $b {break}
    set r0 [lindex $latPos $a0]
    set r1 [lindex $latPos $a1]
    
    # Get the vector between the atoms.
    set d [vecsub $r1 $r0]
    foreach {x y z} $d {break}

    # Map it onto the home cell.
    while {$x < 0.5} {set x [expr $x + 1.0]}
    while {$y < 0.5} {set y [expr $y + 1.0]}
    while {$z < 0.5} {set z [expr $z + 1.0]}
    while {$x >= 0.5} {set x [expr $x - 1.0]}
    while {$y >= 0.5} {set y [expr $y - 1.0]}
    while {$z >= 0.5} {set z [expr $z - 1.0]}
        
    # Map it back into world space.
    set bondVector [vecTransform $basis [list $x $y $z]]

    lappend bondLength [list $a0 $a1 [veclength $bondVector]]
}

# Sort the bond lengths.
foreach zero {0} {set bondLength [lsort -real -index 2 $bondLength]}

# Write the results.
set out [open $outFile w]
foreach b $bondLength {
    puts $out $b
}
close $out
puts "Wrote [llength $bondLength] bonds to $outFile."
exit



