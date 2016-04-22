# Build a coarse-grained SiO-SiN-SiO block.
# to use: vmd -dispdev text -e createSandwich.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Parameters:
set bx 56
set by 64
set bz 24
set dist 4.0
set outerSegname "SIO"
set outerResname "SIO"
set outerName "SIO"
set outerType "PSI"
set outerCharge "  0.000000"
set innerSegname "SIN"
set innerResname "SIN"
set innerName "SIN"
set innerType "CSI"
set innerCharge "  0.000000"
set wrapZ 0
set wrapNone 0

# Output:
set pdb hcp56.pdb
set psf hcp56.psf

set nCells [expr $bx*$by*$bz]
set btop [expr $bz/3]
set bbot [expr 2*$bz/3]
set mass "       72.0000"
set dummy "           0"
set base36 "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"

# Write the psf file.
puts "Writing psf file..."
set out [open $psf w]

# Write the header.
puts $out "PSF"
puts $out ""
puts $out "       1 !NTITLE"
puts $out " REMARKS original generated structure x-plor psf file"

##### ATOMS
puts "Writing atom records..."
puts $out ""
puts $out "[format %8i $nCells] !NATOM"

set n 1

set outerSegCount 0
set outerResid 1
for {set k 0} {$k < $btop} {incr k} {
	for {set j 0} {$j < $by} {incr j} {
		for {set i 0} {$i < $bx} {incr i} {
			set segname [format "%s%s" $outerSegname \
					[string index $base36 $outerSegCount]]
			
			puts -nonewline $out [format "%8i " $n]
			puts -nonewline $out [format "%-4s " $segname]
			puts -nonewline $out [format "%-4i " $outerResid]
			puts -nonewline $out [format "%-4s " $outerResname]
			puts -nonewline $out [format "%-4s " $outerName]
			puts -nonewline $out [format "%-4s " $outerType]
			puts -nonewline $out $outerCharge
			puts -nonewline $out $mass
			puts $out $dummy
			
			incr n
			incr outerResid
			
			# Create a new segment if this one is full.
			if {$outerResid > 9999} {
	    			set outerResid 1
				incr outerSegCount
			}
		}
	}
}

set innerSegCount 0
set innerResid 1
for {set k $btop} {$k < $bbot} {incr k} {
	for {set j 0} {$j < $by} {incr j} {
		for {set i 0} {$i < $bx} {incr i} {
			set segname [format "%s%s" $innerSegname \
					[string index $base36 $innerSegCount]]
			
			puts -nonewline $out [format "%8i " $n]
			puts -nonewline $out [format "%-4s " $segname]
			puts -nonewline $out [format "%-4i " $innerResid]
			puts -nonewline $out [format "%-4s " $innerResname]
			puts -nonewline $out [format "%-4s " $innerName]
			puts -nonewline $out [format "%-4s " $innerType]
			puts -nonewline $out $innerCharge
			puts -nonewline $out $mass
			puts $out $dummy
			
			incr n
			incr innerResid
			
			# Create a new segment if this one is full.
			if {$innerResid > 9999} {
	    			set innerResid 1
				incr innerSegCount
			}
		}
	}
}

for {set k $bbot} {$k < $bz} {incr k} {
	for {set j 0} {$j < $by} {incr j} {
		for {set i 0} {$i < $bx} {incr i} {
			set segname [format "%s%s" $outerSegname \
					[string index $base36 $outerSegCount]]
			
			puts -nonewline $out [format "%8i " $n]
			puts -nonewline $out [format "%-4s " $segname]
			puts -nonewline $out [format "%-4i " $outerResid]
			puts -nonewline $out [format "%-4s " $outerResname]
			puts -nonewline $out [format "%-4s " $outerName]
			puts -nonewline $out [format "%-4s " $outerType]
			puts -nonewline $out $outerCharge
			puts -nonewline $out $mass
			puts $out $dummy
			
			incr n
			incr outerResid
			
			# Create a new segment if this one is full.
			if {$outerResid > 9999} {
	    			set outerResid 1
				incr outerSegCount
			}
		}
	}
}
puts $out ""

# Make the bond list for internal nodes.
puts "Determining bond structure..."
set bond {}
for {set k 0} {$k < $bz} {incr k} {
	set dj [expr 2*($k%2)]
	
	for {set j 0} {$j < $by} {incr j} {
		set di [expr $j%2]
		
		for {set i 0} {$i < $bx} {incr i} {
			# Form a list of the neighbors' 3-indices.
			set neigh {}
			
			# Same row.
			lappend neigh [list [expr $i-1] $j $k]
			lappend neigh [list [expr $i+1] $j $k]
			
			# Top row.
			lappend neigh [list [expr $i-1+$di] [expr $j-1] $k]
			lappend neigh [list [expr $i+$di] [expr $j-1] $k]
			
			# Bottom row.
			lappend neigh [list [expr $i-1+$di] [expr $j+1] $k]
			lappend neigh [list [expr $i+$di] [expr $j+1] $k]
			
			# Outer triangle.
			lappend neigh [list $i $j [expr $k-1]]
			lappend neigh [list [expr $i-1+$di] [expr $j-1+$dj] [expr $k-1]]
			lappend neigh [list [expr $i+$di] [expr $j-1+$dj] [expr $k-1]] 
			
			# Inner triangle.
			lappend neigh [list $i $j [expr $k-1]]
			lappend neigh [list [expr $i-1+$di] [expr $j-1+$dj] [expr $k+1]]
			lappend neigh [list [expr $i+$di] [expr $j-1+$dj] [expr $k+1]]
			
			# Wrap the 3-indices on the periodic boundaries.
			set neighbor {}
			foreach n $neigh {
				foreach {ni nj nk} $n {break}
				if {$wrapNone} {
					if {$ni>=0 && $ni<$bx && $nj>=0 && $nj<$by \
					&& $nk>=0 && $nk<$bz} {
						lappend neighbor [list $ni $nj $nk]
					}
				} elseif {$wrapZ} {
					lappend neighbor [list [expr ($ni+$bx)%$bx] \
					[expr ($nj+$by)%$by] [expr ($nk+$bz)%$bz]]
				} else {
					if {$nk >= 0 && $nk < $bz} {
						lappend neighbor [list [expr ($ni+$bx)%$bx] \
						[expr ($nj+$by)%$by] $nk]
					}
				}
			}
			
			# Map the 3-indices to the atom index.
			set me [expr 1 + $i + $j*$bx + $k*$bx*$by]
			set index {}
			foreach n $neighbor {
				foreach {ni nj nk} $n {break}
				lappend index [expr 1 + $ni + $nj*$bx + $nk*$bx*$by]
			}
			
			# Add the bonds in sorted order.
			foreach n $index {
				lappend bond [lsort -integer [list $me $n]]
			}
		}
	}
}

# Remove redundant bonds.
puts "Removing redundant bonds..."
set totalBonds [llength $bond]
while {1} {
	set bond [lsort -unique $bond]
	break
}
set totalBonds [llength $bond]
puts "Found $totalBonds bonds."

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

# Write the angles.
set totalAngles 0
puts $out ""
puts $out "[format %8i $totalAngles] !NTHETA: angles"
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

##### NON-BONDED
set nNB 0
puts $out ""
puts $out "[format %8i $nNB] !NNB"
puts $out ""

set tmp [expr int($nCells/8)]
set tmp2 [expr $nCells -$tmp*8]
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
puts "The psf file was written successfully."


puts "Generating the pdb file..."
set out [open $pdb w]
set in [open $psf r]
set sx $dist
set sy [expr 0.5*sqrt(3.)*$dist]
set sz [expr sqrt(6.)/3.*$dist]

# Read the atom records from the psf and transfer to the pdb file.
set reading 0
set n 0
set resid0 -1
set serial 1
foreach line [split [read $in] "\n"] {
	if {$reading} {
		# Quit if we have left the atoms section of the psf.
		if {[string match "*!NBOND*" $line]} {break}
		# Check that this is a valid line.
		if {[string length $line] < 33} {continue}
		
		# Find the lattice indices.
		set i [expr $n%$bx]
		set j [expr ($n/$bx)%$by]
		set k [expr ($n/$bx/$by)%$bz]
		set di [expr 0.5*($j%2)]
		set dj [expr 2.*($k%2)/3.]
		
		# Get the atom position.
		set x [expr $sx*($i+$di)]
		set y [expr $sy*($j+$dj)]
		set z [expr $sz*$k]
		
		# Extract the atom data.
		set segname [string range $line 9 12]
		set resid [string trim [string range $line 14 17]]
		set resname [string range $line 19 22]
		set name [string range $line 24 27]
		set type [string range $line 29 32]
		set element "SI"
		
		# Is this a new residue?		
		if {$resid != $resid0} {set serial 1}
						
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
		incr serial
	} elseif {[string match "*!NATOM*" $line]} {
		set reading 1
	}
	
}
close $in
puts "The pdb file was written successfully."

exit





