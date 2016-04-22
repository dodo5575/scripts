# Read the unit cell of a pdb and replicate n1 by n2 by n3 times.
# Use with: vmd -dispdev text -e replicateCrystal.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set nx 4
set ny 4
set nz 1
set basisX {25 0 0}
set basisY {0 25 0}
set basisZ {0 0 55}
# Input:
set name anneal
set psf ${name}_slab.psf
set pdb ${name}_slab.pdb
# Output:
set outPrefix ${name}_surface0

source vector.tcl


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

# Extract the atom records from the record psf.
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
	
	# Skip to the atom records.
	if {[string match "*!NBOND*" $line]} {set reading 1}
    }
    close $in
    return $bondList
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
	set s "[string range $segName0 0 [expr $sl-1]][string index $base36 [expr $index%36]]    "
    } else {
	set s "[string range $segName0 0 2][string index $base36 [expr $index%36]]    "
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

# Make bonds {b0 b1 jx jy jz}, where $jx, $jy, and $jz are image displacements.
proc imageBonds {posVar bondListVar basis} {
    upvar $posVar pos
    upvar $bondListVar bondList
    set basisInv [matInvert $basis]

    set ret {}
    foreach bond $bondList {
	set r0 [lindex $pos [lindex $bond 0]]
	set r1 [lindex $pos [lindex $bond 1]]
	set d [vecSub $r1 $r0]

	set ld [vecTransform $basisInv $d]
	foreach {lx ly lz} $ld {break}
	
	# Find the images that this should connect to.
	set jx 0
	if {$lx < -0.5} {set jx 1}
	if {$lx > 0.5} {set jx -1}
	set jy 0
	if {$ly < -0.5} {set jy 1}
	if {$ly > 0.5} {set jy -1}
	set jz 0
	if {$lz < -0.5} {set jz 1}
	if {$lz > 0.5} {set jz -1}
	
	lappend ret [list [lindex $bond 0] [lindex $bond 1] $jx $jy $jz]
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

# Replicate the system and make the bonds.
proc main {} {
    global psf pdb outPrefix
    global nx ny nz basisX basisY basisZ
    
    set basis [matTranspose [list $basisX $basisY $basisZ]]
    set nImages [expr $nx*$ny*$nz]
    puts "There are $nImages images."
    if {$nImages > 36} {
	puts "Error! Image overflow!"
	exit
    }
    
    set outPsf [open $outPrefix.psf w]
    set outPdb [open $outPrefix.pdb w]
    
    # Write the psf header.
    puts $outPsf "PSF"
    puts $outPsf ""
    puts $outPsf "       7 !NTITLE"
    puts $outPsf " REMARKS original generated structure x-plor psf file"
    puts $outPsf " REMARKS original generated coordinate pdb file"
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
    set nTotal [string range "         [expr $nAtoms]" end-7 end]
    puts $outPsf ""
    puts $outPsf "$nTotal !NATOM"

    # Write the atoms for each image.
    set index 1
    puts "Writing the psf and pdb atom records..."
    for {set j 0} {$j < $nImages} {incr j} {
	puts "Creating image $j"
	set ix [expr $j/$nz/$ny]
	set iy [expr ($j/$nz)%$ny]
	set iz [expr $j%$nz]

	set r0 [vecTransform $basis [list $ix $iy $iz]]
	#set segName [makeSegName "HOG" $j]
	#puts $segName

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

    # Get the bonds.
    set bond [extractPsfBonds $psf]
    set bondImg [imageBonds pos bond $basis]

    # Write the bonds for each image.
    puts "Writing the bonds..."
    set bondLine 0
    set nBonds [expr $nImages*[llength $bondImg]]
    set nBonds1 [string range "         $nBonds" end-7 end]
    puts $outPsf ""
    puts $outPsf "$nBonds1 !NBOND: bonds"
    for {set j 0} {$j < $nImages} {incr j} {
	set ix [expr $j/$nz/$ny]
	set iy [expr ($j/$nz)%$ny]
	set iz [expr $j%$nz]
	puts "Bonding image: $j"

	foreach b $bondImg {
	    foreach {b0 b1 jx jy jz} $b {break}

	    # Find the neighboring image.
	    set kx [expr $ix + $jx]
	    if {$kx < 0} {set kx [expr $kx + $nx]}
	    if {$kx >= $nx} {set kx [expr $kx - $nx]}
	    set ky [expr $iy + $jy]
	    if {$ky < 0} {set ky [expr $ky + $ny]}
	    if {$ky >= $ny} {set ky [expr $ky - $ny]}
	    set kz [expr $iz + $jz]
	    if {$kz < 0} {set kz [expr $kz + $nz]}
	    if {$kz >= $nz} {set kz [expr $kz - $nz]}
	    set k [expr $kz + $ky*$nz + $kx*$nz*$ny]

	    set bond0 [string range "         [expr 1 + $b0 + $j*$n]" end-7 end]
	    set bond1 [string range "         [expr 1 + $b1 + $k*$n]" end-7 end]

	    if {$bondLine < 3} {
		puts -nonewline $outPsf "${bond0}${bond1}"
		incr bondLine
	    } else {
		puts $outPsf "${bond0}${bond1}"
		set bondLine 0
	    }
	}
    }
    
    # Write the rest of the psf crap.
    puts $outPsf ""
    writePsfOther $outPsf $nAtoms
    close $outPsf
    puts "The files $outPrefix.psf and $outPrefix.pdb were written successfully."
}

main
exit


