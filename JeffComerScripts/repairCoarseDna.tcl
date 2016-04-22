# Fix problems introduced by CGBuilder, like incorrect
# and residue numbering.
# Use with: vmd -dispdev text -e repairCoarseDna.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set phoExp PHO
set sugExp SUG
set baseText "(resname ADE GUA and name N1) or (resname CYT THY and name N3)"
#Input:
set atomPsf hairpin.psf
set atomCoords hairpin.pdb
set coarsePdb cg_hairpin.pdb
set database dna.cgc
#Output:
set outPdb cg_hairpin_correct.pdb

proc readDatabase {dataFile} {
    set in [open $dataFile r]
    
    set beadData {}
    set atomData {}
    set nBeadDefs 0
    set n -1

    foreach line [split [read $in] \n] {
	if {[string match "CGBEGIN*" $line]} {
	    set n 0
	    set atom {}
	} elseif {[string match "CGEND*" $line]} {
	    # Add the new bead definition.
	    lappend atomData $atom
	    incr nBeadDefs
	} else {
	    set tok [concat $line]
	    if {[llength $tok] == 3} {
		
		# Read a database line.
		if {$n == 0} {
		    # Extract the bead name.
		    lappend beadData $tok
		    incr n
		} else {
		    # Extract an atom composing the bead.
		    lappend atom $tok
		    incr n
		}
	    }
	}
    }
    close $in

    return [list $beadData $atomData]
}

proc printDatabase {beadData atomData} {
    foreach b $beadData a $atomData {
	puts "($b)"
	foreach atom $a {
	    puts $atom
	}
	puts ""
    }
}

proc centerBeads {coarseMol atomMol beadLine atomLines} {
    set coarseSel [atomselect $coarseMol "name [lindex $beadLine 1]"]
    set res [$coarseSel get resid]
    set o0 [lindex $beadLine 2]

    foreach r $res {
	# Form the selection text for the atoms in this bead.
	set selText ""
	foreach a $atomLines {
	    set name [lindex $a 1]
	    set offset [lindex $a 2]

	    set selText "$selText or (resid [expr $r+$offset-$o0] and name $name)"
	}

	# Remove the initial " or ".
	set selText [string range $selText 4 end]

	# Get the center of mass and shift the bead.
	set sel [atomselect $atomMol $selText]
	if {[$sel num] == 0} {
	    puts "WARNING: Failed to find atoms with $selText."
	} else {
	    set center [measure center $sel weight mass]
	    set bead [atomselect $coarseMol "resid $r and name [lindex $beadLine 1]"]
	    $bead set {x y z} [list $center]
	    $bead delete
	}
	$sel delete
    }
}

# Load the database information.
foreach zero {0} {
    set data [readDatabase $database]
    set beadData [lindex $data 0]
    set atomData [lindex $data 1]
}

# Find the database info refering to phosphates and sugars.
for {set i 0} {$i < [llength $beadData]} {incr i} {
    set name [lindex $beadData $i 1]
    if {[regexp $phoExp $name]} {
	set phoIndex $i
    } elseif {[regexp $sugExp $name]} {
	set sugIndex $i
    }
}

#printDatabase $beadData $atomData

# Load the coarse and atomistic representations.
set coarseMol [mol load pdb $coarsePdb]
set atomMol [mol load psf $atomPsf]
mol addfile $atomCoords

# Reposition the phosphates.
puts "Repositioning phosphate beads..."
centerBeads $coarseMol $atomMol [lindex $beadData $phoIndex] [lindex $atomData $phoIndex]
# Reposition the sugar.
puts "Repositioning sugar beads..."
centerBeads $coarseMol $atomMol [lindex $beadData $sugIndex] [lindex $atomData $sugIndex]
# Reposition the bases.
puts "Repositioning nitrogen bases..."
centerBeads $coarseMol $atomMol {ADE AB 0} {{ADE N1 0}}
centerBeads $coarseMol $atomMol {GUA GB 0} {{GUA N1 0}}
centerBeads $coarseMol $atomMol {CYT CB 0} {{CYT N3 0}}
centerBeads $coarseMol $atomMol {THY TB 0} {{THY N3 0}}

puts "Matching residue IDs to residue names..."
# Get the residue names.
set all [atomselect $coarseMol all]
set sel [atomselect $coarseMol "name \"\[AGCT\]B\""]
foreach zero {0} {
    set res [$sel get resid]
    set resName [$sel get resname]
}
$sel delete
# Enforce correspondence with the residue IDs.
foreach r $res n $resName {
    set sel [atomselect $coarseMol "resid $r"]
    $sel set resname $n
    $sel delete
}

# Remove the extra phosphate.
set last [lindex [lsort -integer [$all get resid]] end]
$all delete
puts "Removing the extra phosphate with residue ID $last..."
set all [atomselect $coarseMol "not resid $last"]

# Write the results.
$all writepdb $outPdb
$all delete
exit



