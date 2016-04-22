# Rename the given system to make it an SiO/SiN/SiO sandwich.
# to use: vmd -dispdev text -e renameAtoms.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set psf0 water_slab1b.psf
set pdb0 water_slab1a.pdb
# Output:
set psf sandwich.psf
set pdb sandwich.pdb

# Parameters:
set mass "       72.0000"
set resname0 "SIO"
set name0 "SIO"
set type0 "PSI"
set charge0 "  0.000000"
set resname1 "SIN"
set name1 "SIN"
set type1 "CSI"
set charge1 "  0.000000"

proc main {} {
	global psf0 pdb0 psf pdb
	global mass
	global resname0 name0 type0 charge0
	global resname1 name1 type1 charge1
	
	set dummy "           0"
	puts "Beginning renameAtoms..."
	
	# Obtain the atom positions.
	puts "Getting the atom positions from $pdb0..."
	set r [extractCoords $pdb0]
	set n [llength $r]
	puts "Found the positions of $n atoms."
	
	set zminmax [getZBounds $r]
	set dz [expr ([lindex $zminmax 1] - [lindex $zminmax 0])/3.]
	set z1 [expr [lindex $zminmax 0] + $dz]
	set z2 [expr $z1 + $dz]
	
	# Rename everything in the psf.
	puts "Generating the psf file using $psf0..."
	set in [open $psf0 r]
	set out [open $psf w]
	set record 0
	set i 0
	foreach line [split [read $in] \n] {
		if {[expr [string length $line] < 2]} {set record 0}
		
		if {$record} {
			set z [lindex $r $i 2]
			set sec [expr $z > $z1 && $z < $z2]
			
			# Extract any needed atom data.
			set segname [string trim [string range $line 9 12]]
			set resid [string trim [string range $line 14 17]]
			
			puts -nonewline $out [format "%8i " [expr $i+1]]
			puts -nonewline $out [format "%-4s " $segname]
			puts -nonewline $out [format "%-4i " $resid]
			puts -nonewline $out [format "%-4s " [set resname${sec}]]
			puts -nonewline $out [format "%-4s " [set name${sec}]]
			puts -nonewline $out [format "%-4s " [set type${sec}]]
			puts -nonewline $out [format "%-4s " [set charge${sec}]]
			puts -nonewline $out $mass
			puts $out $dummy
			
			incr i
		} else {
			# Just write the line if we aren't in the atom records.
			puts $out $line
		}
		
		# Have we reached the atom records?
		if {[string match "*!NATOM*" $line]} {
			puts "Writing atom records..."
			set record 1
		}
	}
	close $in
	close $out
	puts "The psf file $psf was written successfully."
	
	# Rename everything in the pdb.
	puts "Generating the pdb file..."
	makePdbFromPsf $pdb $psf $r
	puts "The pdb file $pdb was written successfully."
}

proc makePdbFromPsf {pdbFile psfFile pos} {
	set out [open $pdbFile w]
	set in [open $psfFile r]

	# Read the atom records from the psf and transfer to the pdb file.
	set n 0
	set reading 0
	foreach line [split [read $in] "\n"] {
		if {$reading} {				
			# Quit if we have left the atoms section of the psf.
			if {[string match "*!NBOND*" $line]} {break}
			# Check that this is a valid line.
			if {[string length $line] < 33} {continue}
		
			# Prepare the coordinates.
			foreach {x y z} [lindex $pos $n] {break}
				
			# Extract the atom data.
			set segname [string range $line 9 12]
			set resid [string trim [string range $line 14 17]]
			set resname [string range $line 19 22]
			set name [string range $line 24 27]
			set type [string range $line 29 32]
			set element "SI"
						
			# Write the pdb line.
			puts -nonewline $out "ATOM  "
			puts -nonewline $out [format "%5i " [expr $n+1]]
			puts -nonewline $out [format " %3s" $name]
			puts -nonewline $out [format "%-4s" $resname]
			puts -nonewline $out " "
			puts -nonewline $out [format "%4s    " $resid]
			puts -nonewline $out [format "%8.3f" $x]
			puts -nonewline $out [format "%8.3f" $y]
			puts -nonewline $out [format "%8.3f" $z]
			puts -nonewline $out "  1.00"
			puts -nonewline $out "  0.00"
			puts -nonewline $out [format "      %-4s" $segname]
			puts -nonewline $out [format "%2s" $element]
			puts $out ""
					
			incr n
		}
		
		# Skip to the atom records.
		if {[string match "*!NATOM*" $line]} {set reading 1}
	}
		
	close $in
	close $out
}

# Returns a list with atom positions.
proc extractCoords {pdbFile} {
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

proc getZBounds {pos} {
	set zmin [lindex $pos 0 2]
	set zmax $zmin
	
	foreach r $pos {
		set z [lindex $r 2]
		if {[expr $z < $zmin]} {set zmin $z}
		if {[expr $z > $zmax]} {set zmax $z}
	}
	
	return [list $zmin $zmax]
}

main
exit





