# Use with: vmd -dispdev text -e tile.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set name ADNA
set segName $name
set periodic 1
set res0 1
set res1 40
# Input:
set psf $name.psf
set pdb $name.pdb
set topFile top_charmm_dna.txt
# Output:
set outPrefix ${name}_charmm

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
	
	# Skip to the atom records.
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
	
	# Skip to the atom records.
	if {[string match "*!NIMPHI*" $line]} {set reading 1}
    }
    close $in
    return $atomList
}

# Construct a psf line from a template line, index, and segname.
proc makePsfLine {template index segName} {
    set newIndex [string range "         $index" end-7 end]
    set newSegName [string range "$segName     " 0 3]
    set tempEnd [string range $template 13 end]
    
    return "${newIndex} ${newSegName}${tempEnd}"
}


# Write the empty angles, dihedrals, impropers, and other stuff.
proc writePsfOther {outPsf nAtoms} {
    puts $outPsf ""
    puts $outPsf "       0 !NTHETA: angles"
    puts $outPsf "\n"
    puts $outPsf "       0 !NPHI: dihedrals"
    puts $outPsf "\n"
    puts $outPsf "       0 !NIMPHI: impropers"
    puts $outPsf ""
    writePsfUnused $outPsf $nAtoms
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

# Write the empty bonds, angles, dihedrals, impropers, and other stuff.
proc writePsfZero {outPsf nAtoms} {
    puts $outPsf ""
    puts $outPsf "       0 !NBOND: bonds"
    puts $outPsf "\n"
    puts $outPsf "       0 !NTHETA: angles"
    puts $outPsf "\n"
    puts $outPsf "       0 !NPHI: dihedrals"
    puts $outPsf "\n"
    puts $outPsf "       0 !NIMPHI: impropers"
    puts $outPsf ""
    writePsfUnused $outPsf $nAtoms
}

proc writeHeader {nRemarks outPsf} {
    puts $outPsf "PSF"
    puts $outPsf ""
    set n [string range "        $nRemarks" end-7 end]
    puts $outPsf "$n !NTITLE"
    return
}

proc writeRemark {remark outPsf} {
    puts $outPsf " REMARKS $remark"
    return
}

proc writeHeaderAtoms {nAtoms outPsf} {
    set nBonds1 [string range "         $nAtoms" end-7 end]
    puts $outPsf "$nBonds1 !NATOM"
    return
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
proc writeItems {itemListVar nRowItems outPsf} {
    upvar $itemListVar itemList

    set rowItem 0
    foreach item $itemList {
	# Insert a newline if necessary.
	if {$rowItem >= $nRowItems} {
	    puts $outPsf ""
	    set rowItem 0
	}

	foreach i $item {
	    # Write the index.
	    set ind [string range "         [expr $i+1]" end-7 end]
	    puts -nonewline $outPsf $ind
	}; # index loop
	incr rowItem
    }; # item loop

    return
}

# Reads the topology into an array topVar(resName,itemType),
# where item type is {bond,angle,dihe,impr}.
proc readTopology {topVar fileName} {
    upvar $topVar top

    # Initialize.
    set res NONE
    set top($res,bond) {}
    set top($res,angle) {}
    set top($res,dihe) {}
    set top($res,impr) {}
    
    set resList {}    
    set in [open $fileName r]
    foreach line [split [read $in] \n] {
	# Check that this is a valid line.
	if {[string length $line] < 2 || [string match "!*" $line]} {continue}
	
	set tok [concat $line]
	if {[string match "RESI *" $line]} {
	    set res [lindex $tok 1]
	    set top($res,bond) {}
	    set top($res,angle) {}
	    set top($res,dihe) {}
	    set top($res,impr) {}
	    lappend resList $res
	} elseif {[string match "BOND *" $line]} {
	    lappend top($res,bond) [lrange $tok 1 2]
	} elseif {[string match "ANGLE *" $line]} {
	    lappend top($res,angle) [lrange $tok 1 3]
	} elseif {[string match "DIHE *" $line]} {
	    lappend top($res,dihe) [lrange $tok 1 4]
	} elseif {[string match "IMPR *" $line]} {
	    lappend top($res,impr) [lrange $tok 1 4]
	} 
    }
    close $in
    return $resList
}

proc linkResidues {topItemVar seg0 res0 seg1 res1} {
    upvar $topItemVar topItem

    set ret {}
    foreach item $topItem {
	set indList {}
	foreach a $item {
	    if {[string match "+*" $a]} {
		set s $seg1
		set r $res1
		set n [string range $a 1 end]
	    } else {
		set s $seg0
		set r $res0
		set n $a
	    }
	    
	    set sel [atomselect top "segname $s and resid $r and name $n"]
	    if {[$sel num] == 0} {
		puts "Error! atom ($s $r $n) does not exist." 
		continue
	    }
	    lappend indList [lindex [$sel get index] 0]
	    $sel delete
	}
	
	lappend ret $indList
    }
    return $ret
}

proc writePsf {fileName remarkVar atomVar bondVar angleVar diheVar imprVar} {
    upvar $remarkVar remark
    upvar $atomVar atom
    upvar $bondVar bond
    upvar $angleVar angle
    upvar $diheVar dihe
    upvar $imprVar impr

    set nRemark [llength $remark]
    set nAtom [llength $atom]
    set nBond [llength $bond]
    set nAngle [llength $angle]
    set nDihe [llength $dihe]
    set nImpr [llength $impr]

    # Open the psf.
    set outPsf [open $fileName w]
    
    # Write the header.
    writeHeader $nRemark $outPsf
    foreach line $remark {
	writeRemark $line $outPsf
    }
    puts $outPsf ""

    # Write the atoms.
    writeHeaderAtoms $nAtom $outPsf
    foreach line $atom {
	puts $outPsf $line
    }
    puts $outPsf ""
    
    # Write the bonds.
    writeHeaderBonds $nBond $outPsf
    writeItems bond 4 $outPsf
    puts $outPsf "\n"
    
    # Write the angles.
    writeHeaderAngles $nAngle $outPsf
    writeItems angle 3 $outPsf
    puts $outPsf "\n"
   
    # Write the dihedrals.
    writeHeaderDihedrals $nDihe $outPsf
    writeItems dihe 2 $outPsf
    puts $outPsf "\n"

    # Write the impropers.
    writeHeaderImpropers $nImpr $outPsf
    writeItems impr 2 $outPsf
    puts $outPsf "\n"

    # Write the unused portions of the psf.
    writePsfUnused $outPsf $nAtom
    close $outPsf
    return 
}

# Replicate the system and make the bonds, angles, dihedrals, impropers.
proc main {} {
    global psf pdb segName res0 res1 topFile outPrefix periodic

    # Extract the topology.
    set resList [readTopology top $topFile]
    puts "Read the topology file $topFile."
    foreach res $resList {
	puts "\nresidue $res:"
	puts "[llength $top($res,bond)] bonds"
	puts "[llength $top($res,angle)] angles"
	puts "[llength $top($res,dihe)] dihedrals"
	puts "[llength $top($res,impr)] impropers"

    }

    # Load the molecule.
    set bond {}
    set angle {}
    set dihe {}
    set impr {}
    mol load psf $psf pdb $pdb
    for {set i $res0} {$i <= $res1} {incr i} {
	set r0 $i
	set r1 [expr $i+1]
	if {$r1 > $res1} {
	    if {$periodic} {
		set r1 $res0
	    } else {
		puts "Warning! Residue $r1 does not exist."
		continue
	    }
	}
	
	set sel [atomselect top "segname $segName and resid $r0"]
	set resName [lindex [$sel get resname] 0]
	if {[$sel num] == 0} {
	    puts "Warning! segname $segName residue $r0 does not exist."
	    continue
	}
	if {[lsearch $resList $resName] < 0} {
	    puts "Warning! no topology for $resName"
	    continue
	}
	

	set bond [concat $bond [linkResidues top($resName,bond) $segName $r0 $segName $r1]]
	set angle [concat $angle [linkResidues top($resName,angle) $segName $r0 $segName $r1]]
	set dihe [concat $dihe [linkResidues top($resName,dihe) $segName $r0 $segName $r1]]
	set impr [concat $impr [linkResidues top($resName,impr) $segName $r0 $segName $r1]]
    }
    
    # Write the psf.
    set remark [list "original generated structure x-plor psf file" "linked from $psf and $pdb by $topFile"]
    set atom [extractPsfAtoms $psf]
    writePsf $outPrefix.psf remark atom bond angle dihe impr
    
    puts "Wrote $outPrefix.psf."
}

main
exit


