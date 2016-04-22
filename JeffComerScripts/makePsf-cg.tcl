# This script will make a psf file for Si3N4.
# use with: vmd -dispdev text -e makePsf-cg.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set pdbFile pore_layer-cg.pdb
# Output:
set fileName [string trimright $pdbFile .pdb]
set boundaryFile ${fileName}.txt
set psfFile ${fileName}.psf
set surfpdb surf.pdb

# "bondDistance" is used to determine whether a bond exists between beads.
set bondDistance 4.8
set numBonds 5
set dummy [format %11g 0]
set type SNCG
set charge 0.
set mass 152.2080
set hexPeriodic 1

# Read the geometry of the system.
set remarkLines {}
set inStream [open $pdbFile r]
foreach line [split [read $inStream] \n] {
    set string0 [string range $line 0 5]

    if { [string match $string0 "REMARK"] } {

	lappend remarkLines $line
	set string1 [string range $line 7 26]

	if { [string match "*unitCell_a*" $string1] } {
	    foreach {tmp1 tmp2 a} $line { break }
	    puts "a = $a"
	}
	if { [string match "*unitCell_b*" $string1] } {
	    foreach {tmp1 tmp2 b} $line { break }
	    puts "b = $b"
	}
	if { [string match "*unitCell_c*" $string1] } {
	    foreach {tmp1 tmp2 c} $line { break }
	    puts "c = $c"
	}
	if { [string match "*basisVec1*" $string1] } {
	    foreach {tmp1 tmp2 x y z} $line { break }
	    set basisVector1 "$x $y $z"
	    puts "basisVec1 = $basisVector1"
	}
	if { [string match "*basisVec2*" $string1] } {
	    foreach {tmp1 tmp2 x y z} $line { break }
	    set basisVector2 "$x $y $z"
	    puts "basisVec2 = $basisVector2"
	}
	if { [string match "*basisVec3*" $string1] } {
	    foreach {tmp1 tmp2 x y z} $line { break }
	    set basisVector3 "$x $y $z"
	    puts "basisVec3 = $basisVector3"
	}
	if { [string match "*nXnYnZ*" $string1] } {
	    foreach {tmp1 tmp2 nX nY nZ} $line { break }
	    puts "nX nY nZ = $nZ $nY $nZ"
	}
    }

}
close $inStream

# Deterimine the lattice vectors.
set vector1 [vecscale $basisVector1 $a]
set vector2 [vecscale $basisVector2 $b]
set vector3 [vecscale $basisVector3 $c]

set pbcVector1 [vecscale $vector1 $nX]
set pbcVector2 [vecscale $vector2 $nY]
set pbcVector3 [vecscale $vector3 $nZ]

puts "Periodic vectors for hexagonal lattice"
set hexVector2 [vecadd [vecscale $vector1 [expr $nX/2]] [vecscale $vector2 [expr $nY/2] ] ]
set hexVector3 [vecadd [vecscale $vector1 [expr $nX/2] ] [vecscale $vector2 [expr -$nY] ] ] 
set hexVector4 [vecadd [vecscale $vector1 [expr $nX]] [vecscale $vector2 [expr -$nY/2] ] ]

puts ""
puts "PERIODIC VECTORS FOR NAMD:"
puts "cellBasisVector1                        $hexVector2"
puts "cellBasisVector2                        [vecinvert $hexVector3]"
puts "cellBasisVector2                        $pbcVector3" 
puts ""

puts "The radius of the hexagon"
set radius [expr 2.*[lindex $hexVector2 0]/3.]

# Write the boundary condition file.
set out [open $boundaryFile w]
puts $out "PERIODIC VECTORS FOR NAMD:"
puts $out "cellBasisVector1                        $hexVector2"
puts $out "cellBasisVector2                        [vecinvert $hexVector3]"
puts $out "cellBasisVector3                        $pbcVector3" 
puts $out ""
puts $out "The radius of the hexagon"
puts $out $radius 
close $out


######## Main Program #########
proc main {} {
global pdbFile
global psfFile
global surfpdb
global bondDistance
global numBonds
global dummy
global type
global charge
global mass
global radius
global hexPeriodic

# Load the pdb.
mol load pdb $pdbFile
set sel [atomselect top all]
set nAtoms [$sel num]
set pos [$sel get {x y z}]
set index [$sel get index] 

# Bond all of the internal atoms.
puts "Bonding internal atoms..."
set bond {}
set bondDistance2 [expr $bondDistance*$bondDistance]
foreach r $pos ind $index {
	# Select neighboring atoms.
	foreach {x y z} $r { break }
	set nearText "($x-x)^2+($y-y)^2+($z-z)^2 < $bondDistance2"
	set near [atomselect top "$nearText and not index $ind"]
	set nearNum [$near num]
	set nearIndex [$near get index]
	
	if {$nearNum > $numBonds} {
		puts "Warning: $nearNum bonds. Reduce bondDistance."
	}
	
	# Add them to the bond list.
	foreach i $nearIndex {lappend bond $ind $i}
}
unset pos
puts "Internal bonds: [expr [llength $bond]/2]"

if {$hexPeriodic} {
# Find the atoms that have fewer than "numBonds" bonds.
puts "Searching for surface atoms..."
set nSurfAtoms 0
foreach ind $index {
	# Find the number of bonds for each atom.
	set n [llength [lsearch -all $bond $ind]]
	set n [expr $n/2]

	# Set the beta value to 0.0 if the atom is on the surface.
	# Set it to 1.0 otherwise.	
	set s [atomselect top "index $ind"]
	if {$n < $numBonds} {
		$s set beta 0.0
		incr nSurfAtoms
	} else {
		$s set beta 1.0
	}
}
puts "Number of surface atoms: $nSurfAtoms"

puts "The system has hexagonal boundary conditions."
# Write a pdb having the surface atoms with beta = 0.0.
$sel writepdb $surfpdb

# Load it up.
mol delete top
mol load pdb $surfpdb

# Determine the centers of the image hexagons.
set sel [atomselect top "beta == 0.0"]
set cen [measure center $sel weight mass]
set pi [expr 4.0*atan(1.0)]
set hexCen {}
set d [expr sqrt(3.)*$radius]
for {set i 0} {$i < 6} {incr i} {
	set theta [expr $pi/6.*(2*$i-1)]
	lappend hexCen [list [expr $d*cos($theta)] [expr $d*sin($theta)] 0.]
}
puts "Periodic image displacements: $hexCen"

# Try to bond surface atoms to the periodic image.
puts "Bonding to the periodic image..."
set surfPos [$sel get {x y z}]
set surfIndex [$sel get index]
foreach hr $hexCen {
	# Shift all of the atoms into this periodic image.
	$sel moveby $hr
	
	foreach ind $surfIndex r $surfPos {
		foreach {x y z} $r {break}
		set nearText "($x-x)^2+($y-y)^2+($z-z)^2 < $bondDistance2"
		
		# Select neighboring atoms.
		set near [atomselect top \
		"beta == 0.0 and $nearText and not index $ind"]
		set nearNum [$near num]
		set nearIndex [$near get index]
		
		# Add them to the bond list.
		foreach i $nearIndex {lappend bond $ind $i}
	}
	
	# Return all atoms to their original position.
	$sel set {x y z} $surfPos
}
unset surfPos
unset surfIndex
}

# Consolidate the bonds into sublists.
puts "Reorganzing bond lists..."
set lBond {}
foreach {b0 b1} $bond {
	lappend lBond [list $b0 $b1]
}

# We should now have all of the bonds twice.
# Find the unique bonds.
puts "Removing redundant bonds..."
set uBond {}
foreach b $lBond {
	set bPerm [list [lindex $b 1] [lindex $b 0]]
	set match [lsearch $uBond $bPerm]
	
	# Add the bond to ubond only if it is unique.
	if {$match == -1} {lappend uBond $b}
}
unset lBond

# Reindex to a 1-based index.
puts "Reindexing..."
set bond {}
foreach b $uBond {
	lappend bond [list [expr [lindex $b 0]+1] [expr [lindex $b 1]+1]]
}
unset uBond
set totalBonds [llength $bond]
#puts $bond
puts "Number of bonds: $totalBonds"

# Find the angles.
puts "Determining the angles..."
set totalBonds1 [expr $totalBonds - 1]
set angle {}
for {set i 0} {$i < $totalBonds1} {incr i} {
	for {set j [expr $i+1]} {$j < $totalBonds} {incr j} {
		foreach {a0 a1} [lindex $bond $i] {break}
		foreach {b0 b1} [lindex $bond $j] {break}
		
		if {$a0 == $b0} {
			lappend angle [list $a1 $a0 $b1]
		} elseif {$a0 == $b1} {
			lappend angle [list $a1 $a0 $b0]
		} elseif {$a1 == $b0} {
			lappend angle [list $a0 $a1 $b1]
		} elseif {$a1 == $b1} {
			lappend angle [list $a0 $a1 $b0]
		}
	}
}
set totalAngles [llength $angle]
puts "Number of angles: $totalAngles"

# Write the psf file.
puts "Writing psf file..."
set out [open $psfFile w]

##### HEADER
puts $out "PSF"
puts $out ""
puts $out "       1 !NTITLE"
puts $out " REMARKS original generated structure x-plor psf file"

##### ATOMS
puts "Writing atom records..."
puts $out ""
puts $out "[format %8i $nAtoms] !NATOM"

set inStream [open $pdbFile r]
set atom 1
foreach line [split [read $inStream] \n] {
    set string0 [string range $line 0 3]

    if { [string match $string0 "ATOM"]  } {
	
	set string1 [string range $line 0 5]
	set string2 [string range $line 6 10]
	set string3 [string range $line 12 15]

	set string4 [string range $line 16 16]

	set string5 [string range $line 17 19]
	set string6 [string range $line 21 21]
	
	set string7 [string range $line 22 25]

	set string8 [string range $line 26 26]

	set string9 [string range $line 30 37]
	set string10 [string range $line 38 45]
	set string11 [string range $line 46 53]

	set string12 [string range $line 54 59]
	set string13 [string range $line 60 65]

	set string14 [string range $line 72 75]

	set string15 [string range $line 76 77]
	set string16 [string range $line 78 79]

	set string3 [string trim $string3]
	puts $out "[format %8i $atom] $string14 [format %-4i $string7] \
	[format %-4s $string5] [format %-4s $string3] [format %-4s $type] \
	$charge $mass $dummy"
	incr atom
    }

}
close $inStream
puts $out ""

##### BONDS
# Write the bonds.
set total [format %8i $totalBonds]
puts $out "$total !NBOND: bonds"
set count 0
foreach b $bond {
	puts -nonewline $out [format "%8i%8i" [lindex $b 0] [lindex $b 1]]
	
	incr count
	if {$count == 4} {
		puts $out ""
		set count 0
	}
}
puts $out ""

##### ANGLES
# Write the angles.
puts $out "[format %8i $totalAngles] !NTHETA: angles"
set count 0
foreach a $angle {
	puts -nonewline $out \
	[format "%8i%8i%8i" [lindex $a 0] [lindex $a 1] [lindex $a 2]]
	
	incr count
	if {$count == 3} {
		puts $out ""
		set count 0
	}
}
puts $out ""

# Write everything else.
##### DIHEDRALS
set nDihedrals 0
puts $out ""
puts $out "[format %8i $nDihedrals] !NPHI: dihedrals"
puts $out ""

##### IMPROPERS
set nImpropers 0
puts $out ""
puts $out "[format %8i $nImpropers] !NIMPHI: impropers"
puts $out ""

##### DONORS
set nDonors 0
puts $out ""
puts $out "[format %8i $nDonors] !NDON: donors"
puts $out ""

##### ACCEPTORS
set nAcceptors 0
puts $out ""
puts $out "[format %8i $nAcceptors] !NACC: acceptors"
puts $out ""

##### NON-BOUNDED
set nNB 0
puts $out ""
puts $out "[format %8i $nNB] !NNB"
puts $out ""

set tmp [expr int($nAtoms/8)]
set tmp2 [expr $nAtoms -$tmp*8]
for {set i 0} {$i <$tmp} {incr i} {
  puts $out "       0       0       0       0       0       0       0       0"
}
set lastString ""
for {set i 0} {$i <$tmp2} {incr i} {
    set lastString "${lastString}       0"
}
puts $out $lastString

####### GROUPS
puts $out ""
puts $out "       1       0 !NGRP"
puts $out "       0       0       0"
puts $out ""
puts $out ""
close $out
}

main
exit



