# Read replicate a psf and pdb nx by ny by nz times.
# The bonds, angles, dihedrals, impropers are reindexed
# appropriately when they cross a replication boundary.
# Use with: vmd -dispdev text -e tile.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Note that unique segment names may not be produced if they 
# have the same first 3 characters!
# Parameters:
set nx 1
set ny 1
set nz 2
set basisX {40 0 0}
set basisY {0 40 0}
set basisZ {0 0 64}
# Input:
set name one_turn_at_pot
set psf ../../${name}.psf
set pdb ../../${name}.pdb
# Output:
set outPrefix ${name}_tiled

source $env(HOME)/scripts/vector.tcl

# Return a list with atom positions.
proc extractPdbCoords {pdbFile} {
    set r {}
    
    # Get the coordinates from the pdb file.
    set in [open $pdbFile r]
    foreach line [split [read $in] \n] {
	if {[string equal [string range $line 0 3] "ATOM"]} {
	    set x [string trim [string range $line 30 37]]
	    set y [string trim [string range $line 38 45]]
	    set z [string trim [string range $line 46 53]]
	    
	    lappend r [list $x $y $z]
	}
    }
    close $in
    return $r
}

# Extract all atom records from a pdb file.
proc extractPdbAtoms {pdbFile} {
    set in [open $pdbFile r]
    
    set pdbLine {}
    foreach line [split [read $in] \n] {
	if {[string equal [string range $line 0 3] "ATOM"]} {
	    lappend pdbLine $line
	}
    }
    close $in	
    
    return $pdbLine
}

# Extract the atom records from the record psf.
proc extractPsfAtoms {recordPsf} {
    set in [open $recordPsf r]
    set reading 0
    set recordList {}
    foreach line [split [read $in] \n] {
	if {$reading} {
	    # Quit if we have left the atoms section of the psf.
	    if {[string match "*!NBOND*" $line]} {break}
	    # Check that this is a valid line.
	    if {[string length $line] < 2} {continue}

	    lappend recordList $line
	}

	# Skip to the atom records.
	if {[string match "*!NATOM*" $line]} {set reading 1}
    }
    close $in
    return $recordList
}

# Extract the bonds from the psf.
# Make them zero-based rather than one-based as in the psf.
proc extractPsfBonds {recordPsf} {
    set in [open $recordPsf r]
    set reading 0
    set bondList {}
    foreach line [split [read $in] \n] {
	if {$reading} {
	    # Quit if we have left the bonds section of the psf.
	    if {[string match "*!NTHETA*" $line]} {break}

	    if {[string equal [string index $line 0] "\#"]} {continue}
	    
	    # Extract the bonds.
	    set tok [concat $line]
	    foreach {b0 b1} $tok {
		lappend bondList [list [expr $b0-1] [expr $b1-1]]
	    }
	}
	
	# Skip to the appropriate section. 
	if {[string match "*!NBOND*" $line]} {set reading 1}
    }
    close $in
    return $bondList
}

# Extract the angles from the psf.
# Make them zero-based rather than one-based as in the psf.
proc extractPsfAngles {recordPsf} {
    set in [open $recordPsf r]
    set reading 0
    set atomList {}
    foreach line [split [read $in] \n] {
	if {$reading} {
	    # Quit if we have left the appropriate section of the psf.
	    if {[string match "*!NPHI*" $line]} {break}

	    if {[string equal [string index $line 0] "\#"]} {continue}
	    
	    # Extract the atom indices.
	    set tok [concat $line]
	    foreach {a0 a1 a2} $tok {
		lappend atomList [list [expr $a0-1] [expr $a1-1] [expr $a2-1]]
	    }
	}
	
	# Skip to the appropriate section.
	if {[string match "*!NTHETA*" $line]} {set reading 1}
    }
    close $in
    return $atomList
}

# Extract the dihedrals from the psf.
# Make them zero-based rather than one-based as in the psf.
proc extractPsfDihedrals {recordPsf} {
    set in [open $recordPsf r]
    set reading 0
    set atomList {}
    foreach line [split [read $in] \n] {
	if {$reading} {
	    # Quit if we have left the appropriate section of the psf.
	    if {[string match "*!NIMPHI*" $line]} {break}

	    if {[string equal [string index $line 0] "\#"]} {continue}
	    
	    # Extract the atom indices.
	    set tok [concat $line]
	    foreach {a0 a1 a2 a3} $tok {
		lappend atomList [list [expr $a0-1] [expr $a1-1] [expr $a2-1] [expr $a3-1]]
	    }
	}
	
	# Skip to the appropriate section. 
	if {[string match "*!NPHI*" $line]} {set reading 1}
    }
    close $in
    return $atomList
}

# Extract the improper dihedrals from the psf.
# Make them zero-based rather than one-based as in the psf.
proc extractPsfImpropers {recordPsf} {
    set in [open $recordPsf r]
    set reading 0
    set atomList {}
    foreach line [split [read $in] \n] {
	if {$reading} {
	    # Quit if we have left the appropriate section of the psf.
	    if {[string match "*!NDON*" $line]} {break}

	    if {[string equal [string index $line 0] "\#"]} {continue}
	    
	    # Extract the atom indices.
	    set tok [concat $line]
	    foreach {a0 a1 a2 a3} $tok {
		lappend atomList [list [expr $a0-1] [expr $a1-1] [expr $a2-1] [expr $a3-1]]
	    }
	}
	
        # Skip to the appropriate section. 
	if {[string match "*!NIMPHI*" $line]} {set reading 1}
    }
    close $in
    return $atomList
}

# Shift a list of vectors by a lattice vector.
proc displaceCell {posVar i1 i2 i3 a1 a2 a3} {
    upvar $posVar pos
    # Compute the new lattice vector.
    set rShift [vecadd [vecscale $i1 $a1] [vecscale $i2 $a2]]
    set rShift [vecadd $rShift [vecscale $i3 $a3]]
    
    set rRep {}
    foreach r $pos {
	lappend rRep [vecadd $r $rShift]
    }
    return $rRep
}

proc getPdbLineSegName {line} {
    return [string range $line 72 75]
}
proc getPsfLineSegName {line} {
    return [string range $line 9 12]
}

proc makeSegName {segName0 index} {
    set base36 "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    set sl [string length [string trim $segName0]]
    if {$sl < 4} {
	set s "[string range $segName0 0 [expr $sl-1]][string index $base36 [expr $index%36]]"
    } else {
	set s "[string range $segName0 0 2][string index $base36 [expr $index%36]]"
    }
    return [string range $s 0 3]
}

# Construct a pdb line from a template line, index, segName, and coordinates.
proc makePdbLine {template index segName r} {
    foreach {x y z} $r {break}
    set record "ATOM  "
    set newIndex [string range "      $index" end-5 end]
    set temp0 [string range $template 12 29]
    set newX [string range [format "       %8.3f" $x] end-7 end]
    set newY [string range [format "       %8.3f" $y] end-7 end]
    set newZ [string range [format "       %8.3f" $z] end-7 end]
    set temp1 [string range $template 54 71]
    set newSegName [string range "$segName     " 0 3]
    set tempEnd [string range $template 76 end]

    # Construct the pdb line.
    return "${record}${newIndex}${temp0}${newX}${newY}${newZ}${temp1}${newSegName}${tempEnd}"
}

# Construct a psf line from a template line, index, and segname.
proc makePsfLine {template index segName} {
    set newIndex [string range "         $index" end-7 end]
    set newSegName [string range "$segName     " 0 3]
    set tempEnd [string range $template 13 end]
    
    return "${newIndex} ${newSegName}${tempEnd}"
}


# Make indices {ai bi ci...} -> {ia 0 0 0 bi bx by bz ci cx cy cz ...}, where .x, .y, and .z are image displacements.
# This works for bonds, angles, dihedrals, impropers.
proc imageNeighbors {posVar atomListVar basis} {
    upvar $posVar pos
    upvar $atomListVar atomList
    set basisInv [matInvert $basis]

    set ret {}
    foreach item $atomList {
	set r0 [lindex $pos [lindex $item 0]]

	# Find the appropriate image for each atom, assuming atom 0
	# is in image 0.
	set newItem {}
	foreach ind $item {
	    set r [lindex $pos $ind]
	    set d [vecSub $r $r0]

	    set ld [vecTransform $basisInv $d]
	    foreach {lx ly lz} $ld {break}

	    # Find the images each atom belongs to.
	    set jx 0
	    if {$lx < -0.5} {set jx 1}
	    if {$lx > 0.5} {set jx -1}
	    set jy 0
	    if {$ly < -0.5} {set jy 1}
	    if {$ly > 0.5} {set jy -1}
	    set jz 0
	    if {$lz < -0.5} {set jz 1}
	    if {$lz > 0.5} {set jz -1}

	    set newItem [concat $newItem [list $ind $jx $jy $jz]]
	}
	lappend ret $newItem
    }
    return $ret
}

# Write the empty angles, dihedrals, impropers, and other stuff.
proc writePsfOther {outPsf nAtoms} {
    puts $outPsf ""
    puts $outPsf "       0 !NTHETA: angles"
    puts $outPsf "\n"
    puts $outPsf "       0 !NPHI: dihedrals"
    puts $outPsf "\n"
    puts $outPsf "       0 !NIMPHI: impropers"
    puts $outPsf "\n"
    puts $outPsf "       0 !NDON: donors"
    puts $outPsf "\n"
    puts $outPsf "       0 !NACC: acceptors"
    puts $outPsf "\n"
    puts $outPsf "       0 !NNB"
    puts $outPsf ""
    set count 0
    for {set i 0} {$i < $nAtoms} {incr i} {
	if {$count == 8} {
	    puts $outPsf ""
	    set count 0
	}
	puts -nonewline $outPsf "       0"
	incr count
    }
    puts $outPsf "\n"
    puts $outPsf "       1       0 !NGRP"
    puts $outPsf "       0       0       0"
    puts $outPsf ""
}

# Write the unused portions of the psf.
proc writePsfUnused {outPsf nAtoms} {
    puts $outPsf ""
    puts $outPsf "       0 !NDON: donors"
    puts $outPsf "\n"
    puts $outPsf "       0 !NACC: acceptors"
    puts $outPsf "\n"
    puts $outPsf "       0 !NNB"
    puts $outPsf ""
    set count 0
    for {set i 0} {$i < $nAtoms} {incr i} {
	if {$count >= 8} {
	    puts $outPsf ""
	    set count 0
	}
	puts -nonewline $outPsf "       0"
	incr count
    }
    puts $outPsf "\n\n"
    puts $outPsf "       1       0 !NGRP"
    puts $outPsf "       0       0       0"
    puts $outPsf ""
}

proc writeHeaderBonds {nBonds outPsf} {
    set nBonds1 [string range "         $nBonds" end-7 end]
    puts $outPsf "$nBonds1 !NBOND: bonds"
    return
}

proc writeHeaderAngles {nAngles outPsf} {
    set nAngles1 [string range "         $nAngles" end-7 end]
    puts $outPsf "$nAngles1 !NTHETA: angles"
    return
}

proc writeHeaderDihedrals {nDihedrals outPsf} {
    set nDihedrals1 [string range "         $nDihedrals" end-7 end]
    puts $outPsf "$nDihedrals1 !NPHI: dihedrals"
    return
}

proc writeHeaderImpropers {nImpropers outPsf} {
    set nImpropers1 [string range "         $nImpropers" end-7 end]
    puts $outPsf "$nImpropers1 !NIMPHI: impropers"
    return
}

# Write a set of items (bonds, angles, dihedrals, impropers)
# with $nRowItems per row.
# $itemListVar has the format {ai ax ay az b bx...} where
# $ai is the atom index and {$ax $ay $az} are the image indices. 
proc writeImageItems {itemListVar imgSize nAtomsPerImage nRowItems outPsf} {
    upvar $itemListVar itemList
    foreach {nx ny nz} $imgSize {break}
    set nImages [expr $nx*$ny*$nz]
    set n $nAtomsPerImage
    
    set rowItem 0
    for {set j 0} {$j < $nImages} {incr j} {
	set ix [expr $j/$nz/$ny]
	set iy [expr ($j/$nz)%$ny]
	set iz [expr $j%$nz]

	foreach item $itemList {
	    # Insert a newline if necessary.
	    if {$rowItem >= $nRowItems} {
		puts $outPsf ""
		set rowItem 0
	    }

	    foreach {ai ax ay az} $item {
		# Find the neighboring image.
		set kx [expr $ix + $ax]
		if {$kx < 0} {set kx [expr $kx + $nx]}
		if {$kx >= $nx} {set kx [expr $kx - $nx]}
		set ky [expr $iy + $ay]
		if {$ky < 0} {set ky [expr $ky + $ny]}
		if {$ky >= $ny} {set ky [expr $ky - $ny]}
		set kz [expr $iz + $az]
		if {$kz < 0} {set kz [expr $kz + $nz]}
		if {$kz >= $nz} {set kz [expr $kz - $nz]}
		set k [expr $kz + $ky*$nz + $kx*$nz*$ny]


		# Write the index.
		set ind [expr 1 + $ai + $k*$n]
		set ind1 [string range "         $ind" end-7 end]
		puts -nonewline $outPsf $ind1
	    }; # index loop

	    incr rowItem
	}; # item loop
    }; # image loop
    return [expr $nImages*[llength $itemList]]
}

# Replicate the system and make the bonds, angles, dihedrals, impropers.
proc main {} {
    global psf pdb outPrefix
    global nx ny nz basisX basisY basisZ
    
    set imageSize [list $nx $ny $nz]
    set basis [matTranspose [list $basisX $basisY $basisZ]]
    set nImages [expr $nx*$ny*$nz]
    puts "There are $nImages images."
    if {$nImages > 36} {
	puts "Error! Too many images. Cannot produce unique segment names."
	exit
    }
    
    set outPsf [open $outPrefix.psf w]
    set outPdb [open $outPrefix.pdb w]
    
    # Write the psf header.
    puts $outPsf "PSF"
    puts $outPsf ""
    puts $outPsf "       6 !NTITLE"
    puts $outPsf " REMARKS original generated structure x-plor psf file"
    puts $outPsf " REMARKS replicated from $psf and $pdb"
    puts $outPsf " REMARKS basisVector1 $basisX" 
    puts $outPsf " REMARKS basisVector2 $basisY"
    puts $outPsf " REMARKS basisVector3 $basisZ" 
    puts $outPsf " REMARKS replicationCount $nx $ny $nz" 
    
    # Write the pdb header.
    puts $outPdb "REMARK original generated coordinate pdb file"
    puts $outPdb "REMARK replicated from $psf and $pdb"
    puts $outPdb "REMARK basisVector1 $basisX" 
    puts $outPdb "REMARK basisVector2 $basisY"
    puts $outPdb "REMARK basisVector3 $basisZ" 
    puts $outPdb "REMARK replicationCount $nx $ny $nz" 

    # Get the psf and pdb lines.
    set psfAtoms [extractPsfAtoms $psf]
    set pdbAtoms [extractPdbAtoms $pdb]
    puts "Read [llength $psfAtoms] atom records from $psf."
    puts "Read [llength $pdbAtoms] atom records from $pdb."

    #puts [lrange $psfAtoms 0 2]
    #puts "..."
    #puts [lrange $psfAtoms end-2 end]

    # Get the coordinates.
    set pos [extractPdbCoords $pdb]
    set n [llength $pos]
    puts "Read $n positions from $pdb."

    # Start the atom records.
    set nAtoms [expr $nImages*$n]
    set nAtoms1 [string range "         [expr $nAtoms]" end-7 end]
    puts $outPsf ""
    puts $outPsf "$nAtoms1 !NATOM"

    # Write the atoms for each image.
    set index 1
    puts "Writing the psf and pdb atom records..."
    for {set j 0} {$j < $nImages} {incr j} {
	puts "Creating image $j"
	set ix [expr $j/$nz/$ny]
	set iy [expr ($j/$nz)%$ny]
	set iz [expr $j%$nz]

	set r0 [vecTransform $basis [list $ix $iy $iz]]
	
	foreach psfLine $psfAtoms pdbLine $pdbAtoms r $pos {
	    set oldSegName [getPdbLineSegName $pdbLine]
	    set segName [makeSegName $oldSegName $j]
	    set r1 [vecAdd $r $r0]

	    puts $outPsf [makePsfLine $psfLine $index $segName]
	    puts $outPdb [makePdbLine $pdbLine $index $segName $r1] 
	    incr index
	}
    }
    close $outPdb
    puts "`$outPrefix.pdb' written successfully."
 
    # Write the bonds for each image.
    set item [extractPsfBonds $psf]
    set itemList [imageNeighbors pos item $basis]
    puts $outPsf ""
    set nTotal [expr $nImages*[llength $itemList]]
    writeHeaderBonds $nTotal $outPsf
    puts "Writing $nTotal bonds..."
    writeImageItems itemList $imageSize $n 4 $outPsf
    puts $outPsf ""
    
    # Write the angles for each image.
    set item [extractPsfAngles $psf]
    set itemList [imageNeighbors pos item $basis]
    puts $outPsf ""
    set nTotal [expr $nImages*[llength $itemList]]
    writeHeaderAngles $nTotal $outPsf
    puts "Writing $nTotal angles..."
    writeImageItems itemList $imageSize $n 3 $outPsf
    puts $outPsf ""
   
    # Write the dihedrals for each image.
    set item [extractPsfDihedrals $psf]
    set itemList [imageNeighbors pos item $basis]
    puts $outPsf ""
    set nTotal [expr $nImages*[llength $itemList]]
    writeHeaderDihedrals $nTotal $outPsf
    puts "Writing $nTotal dihedrals..."
    writeImageItems itemList $imageSize $n 2 $outPsf
    puts $outPsf ""

    # Write the impropers for each image.
    set item [extractPsfImpropers $psf]
    set itemList [imageNeighbors pos item $basis]
    puts $outPsf ""
    set nTotal [expr $nImages*[llength $itemList]]
    writeHeaderImpropers $nTotal $outPsf
    puts "Writing $nTotal impropers..."
    writeImageItems itemList $imageSize $n 2 $outPsf
    puts $outPsf ""

    # Write the unused portions of the psf.
    writePsfUnused $outPsf $nAtoms
    close $outPsf
    puts "`$outPrefix.psf' written successfully."
}

main
exit


