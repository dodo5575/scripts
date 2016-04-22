# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set cutoff 20.0
#Input:
set psf pore_grisha.psf
set pdb pore_grisha.pdb
#Output:
set finalPsf pore_grisha_surf.psf
set finalPdb pore_grisha_surf.pdb

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

# Find bonds apparently longer than the cutoff.
set surfIndex {}
foreach b $bond {
    set d [measure bond [list [lindex $b 0] [lindex $b 1]]]
    if {$d > $cutoff} {
	lappend surfIndex [lindex $b 0]
	lappend surfIndex [lindex $b 1]
    }
}
foreach zero {0} {set surfIndex [lsort -unique -integer $surfIndex]}
puts "\n[llength $surfIndex] surface atoms were found."

# Obtain the {segid resid name} for the selection.
set sel [atomselect top "not (index $surfIndex)"]
puts "[$sel num] interior atoms."
foreach zero {0} {set atomList [lsort -unique [$sel get {segname resid name}]]}
set nAtoms [$sel num]
set nResidues [llength $atomList]
$sel delete

package require psfgen 1.3
resetpsf

readpsf $psf
coordpdb $pdb

# Delete the selection.
foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
}

writepsf $finalPsf
writepdb $finalPdb
puts ""
puts "$nAtoms atoms were deleted."
puts "$nResidues residues were deleted."
mol delete top
exit



