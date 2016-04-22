# Read replicate a psf and pdb nx by ny by nz times.
# The bonds, angles, dihedrals, impropers are reindexed
# appropriately when the cross a replication boundary.
# Use with: vmd -dispdev text -e tile.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>


# Input:
set templateName cont_dna_charmm
set templatePsf ${templateName}.psf
set templatePdb ${templateName}.pdb
set templateRes0 {BDNA 80}
set templateRes1 {BDNA 1}
# Output:
set topFile test_charmm_ade.txt

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

proc coLinks {itemList indList0 indList1} {
    set ret {}
    
    foreach item $itemList {
	# Check that the item contains only indices indList0 union indList1.
	set good 1
	foreach i $item {
	    if {[lsearch $indList0 $i] < 0 && [lsearch $indList1 $i] < 0} {
		set good 0
		break
	    }
	}

	# Check that the list contains at least one index from indList0.
	if {$good} {
	    set good0 0

	    foreach i $item {
		if {[lsearch $indList0 $i] >= 0} {
		    set good0 1
		    break
		}
	    }
	    if {$good0} {
		lappend ret $item
	    }
	}
    }
    return $ret
}

proc internalLinks {itemList indList} {
    set ret {}

    foreach item $itemList {
	set good 1

	foreach i $item {
	    if {[lsearch $indList $i] < 0} {
		set good 0
		break
	    }
	}

	if {$good} {
	    lappend ret $item
	}
    }
    return $ret
}

proc externalLinks {itemList indList0 indList1} {
    set ret {}

    foreach item $itemList {
	set good 1

	# Check that the item contains only indices indList0 union indList1.
	foreach i $item {
	    if {[lsearch $indList0 $i] < 0 && [lsearch $indList1 $i] < 0} {
		set good 0
		break
	    }
	}
	
	# Ignore items that contain only internal links.
	if {$good} {
	    set good0 0
	    set good1 0
	    
	    foreach i $item {
		if {[lsearch $indList0 $i] >= 0} {
		    set good0 1
		}
		if {[lsearch $indList1 $i] >= 0} {
		    set good1 1
		}
	    }
	    
	    if {$good0 && $good1} {
		lappend ret $item
	    }
	}
    }
    return $ret
}

proc topologyLine {item indList0 nameList0 indList1 nameList1} {   
    set line ""
    foreach ind $item {
	set j [lsearch $indList0 $ind]
	if {$j >= 0} {
	    set line "$line [lindex $nameList0 $j]"
	} else {
	    set j [lsearch $indList1 $ind]
	    set line "$line +[lindex $nameList1 $j]"
	}
    }

    return $line
}

# Replicate the system and make the bonds, angles, dihedrals, impropers.
proc writeTopology {templatePsf templatePdb templateRes0 templateRes1 topFile} {
    set res0 [lindex $templateRes0 1]
    set res1 [lindex $templateRes1 1]

    # Get the indices of the template's residues.
    mol load psf $templatePsf pdb $templatePdb
    set sel [atomselect top "segname [lindex $templateRes0 0] and resid $res0"]
    set indList0 [$sel get index]
    set nameList0 [$sel get name]
    set resName [lindex [$sel get resname] 0]
    $sel delete
    
    set sel [atomselect top "segname [lindex $templateRes1 0] and resid $res1"]
    set indList1 [$sel get index]
    set nameList1 [$sel get name]
    $sel delete
    mol delete top

    # Write the topology file.
    set outTop [open $topFile w]
    puts $outTop "RESI $resName"
    set topBonds {}
    set topAngles {}
    set topDihedrals {}
    set topImpropers {}
    
    set itemList [coLinks [extractPsfBonds $templatePsf] $indList0 $indList1]
    foreach item $itemList {
	lappend topBonds [topologyLine $item $indList0 $nameList0 $indList1 $nameList1] 
    }
    puts "Found [llength $topBonds] bonds."

    set itemList [coLinks [extractPsfAngles $templatePsf] $indList0 $indList1]
    foreach item $itemList {
	lappend topAngles [topologyLine $item $indList0 $nameList0 $indList1 $nameList1] 
    }
    puts "Found [llength $topAngles] angles."

    set itemList [coLinks [extractPsfDihedrals $templatePsf] $indList0 $indList1]
    foreach item $itemList {
	lappend topDihedrals [topologyLine $item $indList0 $nameList0 $indList1 $nameList1] 
    }
    puts "Found [llength $topDihedrals] dihedrals."

    set itemList [coLinks [extractPsfImpropers $templatePsf] $indList0 $indList1]
    foreach item $itemList {
	lappend topImpropers [topologyLine $item $indList0 $nameList0 $indList1 $nameList1] 
    }
    puts "Found [llength $topImpropers] impropers."
    
    foreach item [lsort $topBonds] {
	puts $outTop "BOND $item"
    }
    foreach item [lsort $topAngles] {
	puts $outTop "ANGLE $item"
    }
    foreach item [lsort $topDihedrals] {
	puts $outTop "DIHE $item"
    }
    foreach item [lsort $topImpropers] {
	puts $outTop "IMPR $item"
    }

    return 0
}

writeTopology $templatePsf $templatePdb $templateRes0 $templateRes1 $topFile
exit


