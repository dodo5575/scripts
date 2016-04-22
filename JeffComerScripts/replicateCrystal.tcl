# Read the unit cell of a pdb and replicate n1 by n2 by n3 times.
# Use with: vmd -dispdev text -e replicateCrystal.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set unitCellPdb unit_cell_alpha.pdb
# Output:
set outPdb membrane.pdb
# Parameters:
# Choose n1 and n2 even if you wish to use cutHexagon.tcl.
set n1 6
set n2 6
set n3 6
set l1 7.595
set l2 7.595
set l3 2.902
set basisVector1 [list 1.0 0.0 0.0]
set basisVector2 [list 0.5 [expr sqrt(3.)/2.] 0.0]
set basisVector3 [list 0.0 0.0 1.0]

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
proc extractPdbRecords {pdbFile} {
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

# Shift a list of vectors by a lattice vector.
proc displaceCell {rUnitName i1 i2 i3 a1 a2 a3} {
	upvar $rUnitName rUnit
	# Compute the new lattice vector.
	set rShift [vecadd [vecscale $i1 $a1] [vecscale $i2 $a2]]
	set rShift [vecadd $rShift [vecscale $i3 $a3]]
	
	set rRep {}
	foreach r $rUnit {
		lappend rRep [vecadd $r $rShift]
	}
	return $rRep
}

# Construct a pdb line from a template line, index, resId, and coordinates.
proc makePdbLine {template index resId r} {
	foreach {x y z} $r {break}
	set record "ATOM  "
	set si [string range [format "     %5i " $index] end-5 end]
	set temp0 [string range $template 12 21]
	set resId [string range "    $resId"  end-3 end]
	set temp1 [string range $template  26 29]
	set sx [string range [format "       %8.3f" $x] end-7 end]
	set sy [string range [format "       %8.3f" $y] end-7 end]
	set sz [string range [format "       %8.3f" $z] end-7 end]
	set tempEnd [string range $template 54 end]

	# Construct the pdb line.
	return "${record}${si}${temp0}${resId}${temp1}${sx}${sy}${sz}${tempEnd}"
}

# Build the crystal.
proc main {} {
	global unitCellPdb outPdb
	global n1 n2 n3 l1 l2 l3 basisVector1 basisVector2 basisVector3
		
	set out [open $outPdb w]
	puts $out "REMARK Unit cell dimensions:"
	puts $out "REMARK a1 $l1" 
	puts $out "REMARK a2 $l2" 
	puts $out "REMARK a3 $l3" 
	puts $out "REMARK Basis vectors:"
	puts $out "REMARK basisVector1 $basisVector1" 
	puts $out "REMARK basisVector2 $basisVector2" 
	puts $out "REMARK basisVector3 $basisVector3" 
	puts $out "REMARK replicationCount $n1 $n2 $n3" 
	
	set a1 [vecscale $l1 $basisVector1]
	set a2 [vecscale $l2 $basisVector2]
	set a3 [vecscale $l3 $basisVector3]
	
	set rUnit [extractPdbCoords $unitCellPdb]
	set pdbLine [extractPdbRecords $unitCellPdb]
	puts "\nReplicating unit $unitCellPdb cell $n1 by $n2 by $n3..."
		
	# Replicate the unit cell.
	set atom 1
	set resId 1
	for {set k 0} {$k < $n3} {incr k} {
    		for {set j 0} {$j < $n2} {incr j} {
			for {set i 0} {$i < $n1} {incr i} {
				set rRep [displaceCell rUnit $i $j $k $a1 $a2 $a3]
				
				# Write each atom.
				foreach r $rRep l $pdbLine {
					puts $out [makePdbLine $l $atom $resId $r]
					incr atom
				}
				incr resId
				
				if {$resId > 9999} {
					puts "Warning! Residue overflow."
					set resId 1
				}								
			}
		}
	}
	puts $out "END"
	close $out
	
	puts "The file $outPdb was written successfully."
}

main



