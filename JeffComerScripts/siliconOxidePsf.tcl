# Make a psf file for SiO2.
# Use with: vmd -dispdev text -e siliconOxidePsf.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set fileNamePrefix silica_pore
# Input:
set pdbFile ${fileNamePrefix}1.pdb
set boundaryFile ${fileNamePrefix}1.bound
# Output:
set psfFile ${fileNamePrefix}2.psf
set surfPdb surf.pdb
# Parameters:
# Should angles be calculated in addition to bonds?
set findAngles 0
set hexPeriodic 0
set squarePeriodic 1
set zPeriodic 0
# "bondDistance" is used to determine whether a bond exists between atoms.
set bondDistance 2.0
# Si parameters
set nameSi "SI.*"
set massSi 28.085500
set chargeSi 1.0
set typePrefixSi SI_
set numBondsSi 4
# N parameters
set nameN "O.*"
set massN 15.9994
# chargeN is determined by neutrality.
set chargeN -0.5
set typePrefixN O_
set numBondsN 2

proc main {} {
	global pdbFile boundaryFile psfFile surfPdb
	global findAngles hexPeriodic squarePeriodic zPeriodic
	global bondDistance
	global nameSi massSi chargeSi typePrefixSi numBondsSi
	global nameN massN chargeN typePrefixN numBondsN
	
	set selTextSi "name \"${nameSi}\""
	set selTextN "name \"${nameN}\""
	
	# Load the pdb.
	mol load pdb $pdbFile
	set nAtoms [molinfo top get numatoms]
	
	# Get the number of nitrogen and silicon atoms.
	set silicon [atomselect top $selTextSi]
	set numSilicon [$silicon num]
	$silicon delete
	set nitrogen [atomselect top $selTextN]
	set numNitrogen [$nitrogen num]
	$nitrogen delete
	
	# Determine the nitrogen charge.
	set chargeN0 $chargeN
	set chargeN [expr -$chargeSi*$numSilicon/$numNitrogen]
	puts "\nSet charge on $nameN to $chargeN to enforce neutrality."
	puts "Relative error: [expr ($chargeN-$chargeN0)/$chargeN0]\n"
	
	# Find the internal bonds.
	puts "Bonding internal atoms..."
	set bond [bondAtoms all $bondDistance]
	puts "Internal bonds: [expr [llength $bond]/4]"
	
	# Bond to periodic images.	
	if {$hexPeriodic || $squarePeriodic || $zPeriodic} {
		# Create the surface atom pdb.
		set all [atomselect top all]
		$all set beta 1.0
		puts "Searching for surface atoms..."
		set nSurfSi [markSurface $bond $selTextSi $numBondsSi]
		set nSurfN [markSurface $bond $selTextN $numBondsN]
		puts "Number of surface silicons: $nSurfSi"
		puts "Number of surface nitrogens: $nSurfN"
		$all writepdb $surfPdb
		$all delete
				
		# Load it up.
		mol delete top
		mol load pdb $surfPdb
		
		if {$hexPeriodic} {
			puts "The system has hexagonal periodic boundary conditions."
			set radius [readRadius $boundaryFile]
			puts "Hexagon radius: $radius"
			
			# Determine the centers of the image hexagons.
			set pi [expr 4.0*atan(1.0)]
			set hexCen {}
			set d [expr sqrt(3.)*$radius]
			for {set i 0} {$i < 6} {incr i} {
				set theta [expr $pi/6.*(2*$i-1)]
				lappend hexCen [list [expr $d*cos($theta)] \
				[expr $d*sin($theta)] 0.]
			}
			puts "Periodic image displacements: $hexCen"
		
			# Find the bonds on the periodic boundaries.
			puts "Bonding to the periodic image..."
			foreach r $hexCen {
				set bond [concat $bond \
				[bondPeriodic all $bondDistance $r]]
			}
		}
		
		if {$squarePeriodic} {
			puts "The system has square periodic boundary conditions."
			set radius [readRadius $boundaryFile]
			puts "Square radius: $radius"
			
			# Determine the centers of the image squares.
			set pi [expr 4.0*atan(1.0)]
			set squareCen {}
			set d [expr 2*$radius]
			for {set i 0} {$i < 4} {incr i} {
				set theta [expr $pi/2.*$i]
				lappend squareCen [list [expr $d*cos($theta)] \
				[expr $d*sin($theta)] 0.]
			}
			puts "Periodic image displacements: $squareCen"
		
			# Find the bonds on the periodic boundaries.
			puts "Bonding to the periodic image..."
			foreach r $squareCen {
				set bond [concat $bond \
				[bondPeriodic all $bondDistance $r]]
			}
		}
			
		if {$zPeriodic} {
			puts "The system is periodic along the z-axis."
			set lz [readLz $boundaryFile]
			puts "Period in z: $lz"
			set zCen [list [list 0 0 -${lz}] [list 0 0 $lz]]
		
			# Find the bonds on the periodic boundaries.
			puts "Bonding to the periodic image..."
			foreach r $zCen {
				set bond [concat $bond \
				[bondPeriodic all $bondDistance $r]]
			}
		}
	
	} 
	mol delete top
	
	puts "Counting bonds on each atom..."
	countBonds count $bond $nAtoms
	puts "Reorganizing bond lists..."
	set bond [reorganizeBonds $bond]
	puts "Removing redundancy..."
	set bond [removeRedundantBonds $bond]
	set totalBonds [llength $bond]
	puts "Number of bonds: $totalBonds"
	
	set angle {}
	if {$findAngles} {
		puts "Determining the angles..."
		set angle [findAngles $bond]
		set totalAngles [llength $angle]
		puts "Number of angles: $totalAngles"
	}
	
	puts "Writing psf file..."
	manifestPsf $psfFile $pdbFile $nAtoms bond angle count
	puts "The file $psfFile was written successfully."
}

# Find bonds between internal atoms and return them.
proc bondAtoms {selText bondDistance} {
	set sel [atomselect top $selText]
	set pos [$sel get {x y z}]
	set index [$sel get index]
	$sel delete

	set bondDistance2 [expr $bondDistance*$bondDistance]
	set bond {}
	foreach r $pos ind $index {
		# Select neighboring atoms.
		foreach {x y z} $r { break }
		set nearText "($x-x)^2+($y-y)^2+($z-z)^2 < $bondDistance2"
		set near [atomselect top \
		"$selText and $nearText and not index $ind"]
		set nearNum [$near num]
		set nearIndex [$near get index]
		$near delete
	
		# Add them to the bond list.
		foreach i $nearIndex {lappend bond $ind $i}
	}
	return $bond
}

# Get the radius from the boundary file.
proc readRadius {boundaryFile} {
	set in [open $boundaryFile r]
	foreach line [split [read $in] \n] {
    		if {[string match "radius *" $line]} {
			set radius [lindex $line 1]
			break
		}
	}
	close $in
	return $radius
}

# Get the cellBasisVector3_z from the boundary file.
proc readLz {boundaryFile} {
	set in [open $boundaryFile r]
	foreach line [split [read $in] \n] {
    		if {[string match "cellBasisVector3 *" $line]} {
			set lz [lindex $line 3]
			break
		}
	}
	close $in
	return $lz
}

# Try to bond surface atoms to the periodic image.
proc bondPeriodic {selText bondDistance periodicDisp} {
	set selText "$selText and beta == 0.0"
	set sel [atomselect top $selText]
	set pos [$sel get {x y z}]
	set index [$sel get index]
		
	# Shift all of the atoms into this periodic image.
	$sel moveby $periodicDisp
		
	set bondDistance2 [expr $bondDistance*$bondDistance]
	set bond {}
	foreach r $pos ind $index {
		# Select neighboring atoms.
		foreach {x y z} $r { break }
		set nearText "($x-x)^2+($y-y)^2+($z-z)^2 < $bondDistance2"
		set near [atomselect top \
		"$selText and $nearText and not index $ind"]
		set nearNum [$near num]
		set nearIndex [$near get index]
		$near delete
	
		# Add them to the bond list.
		foreach i $nearIndex {lappend bond $ind $i}
	}

	# Return all atoms to their original position.
	$sel set {x y z} $pos
	$sel delete
		
	return $bond
}


# Find the atoms that have fewer than "numBonds" bonds.
# Mark surface atoms by beta = 0.0.
# Warning! The bond list is assumed to be flat and redundant.
proc markSurface {bond selText numBonds} {
	set sel [atomselect top $selText]
	set index [$sel get index]
	set nSurfAtoms 0
	
	foreach i $index {
		# Find the number of bonds for each atom.
		set n [llength [lsearch -all $bond $i]]
		# Assume each bond is in the list twice.
		set n [expr $n/2]
		
		# Set the beta value to 0.0 if the atom is on the surface.
		if {$n < $numBonds} {
			set s [atomselect top "index $i"]
			$s set beta 0.0
			incr nSurfAtoms
			$s delete	
		}
	}
	$sel delete
	
	return $nSurfAtoms
}

# Count the number of bonds on each atom and return an array (zero-based).
# The result is placed in a variable name countVar.
# Warning! The bond list is assumed to be flat and redundant.
proc countBonds {countVar bond nAtoms} {
	upvar $countVar count
	
	set num {}
	for {set i 0} {$i < $nAtoms} {incr i} {
		set n [llength [lsearch -all $bond $i]]
		set n [expr $n/2]
		lappend num $i $n
	}
	
	array set count $num
}

# Put the bonds into sublists.
# Reindex to a 1-based index.
proc reorganizeBonds {bond} {
	set ret {}
	foreach {b0 b1} $bond {
		incr b0
		incr b1
		lappend ret [list $b0 $b1]
	}
	return $ret
}

# We should now have all of the bonds twice.
# Find the unique bonds.
proc removeRedundantBonds {bond} {
	set ret {}
	foreach b $bond {
		set bPerm [list [lindex $b 1] [lindex $b 0]]
		set match [lsearch $ret $bPerm]
	
		# Add the bond to "ret" only if it is unique.
		if {$match == -1} {lappend ret $b}
	}
	return $ret
}

# Find the angles.
proc findAngles {bond} {
	set totalBonds [llength $bond]
	set totalBonds1 [expr $totalBonds - 1]

	# Find bonds that share atoms.
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
	return $angle
}

# Write the psf file.
proc manifestPsf {psfFile pdbFile nAtoms bondVar angleVar countVar} {
	global nameSi massSi chargeSi typePrefixSi numBondsSi
	global nameN massN chargeN typePrefixN numBondsN
	
	# Import the big pass-by-reference stuff.
	upvar $bondVar bond
	upvar $angleVar angle
	upvar $countVar count
		 
	set dummy "          0"
	set totalBonds [llength $bond]
	set totalAngles [llength $angle]
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
	
	# Open the pdb to extract the atom records.
	set inStream [open $pdbFile r]
	set atom 1
	foreach line [split [read $inStream] \n] {
    		set string0 [string range $line 0 3]
		if {![string match $string0 "ATOM"]} {continue}
		
		# Extract each pdb field.
		set record [string range $line 0 5]
		set serial [string range $line 6 10]
		set name [string range $line 12 15]
		set altLoc [string range $line 16 16]
		set resName [string range $line 17 19]
		set chainId [string range $line 21 21]
		set resId [string range $line 22 25]
		set iCode [string range $line 26 26]
		set x [string range $line 30 37]
		set y [string range $line 38 45]
		set z [string range $line 46 53]
		set occupancy [string range $line 54 59]
		set beta [string range $line 60 65]
		set segName [string range $line 72 75]
		set element [string range $line 76 77]
		set charge [string range $line 78 79]
		
		# Determine the type names.
		set numBonds $count([expr $atom-1])
		set typeSi ${typePrefixSi}${numBonds}
		set typeN ${typePrefixN}${numBonds}
		
		# Write the atom record.	
		puts -nonewline $out [format "%8i " $atom]
		puts -nonewline $out [format "%-4s " $segName]
		puts -nonewline $out [format "%-4i " $resId]
		puts -nonewline $out [format "%-3s " $resName]
		puts -nonewline $out [format "%-4s  " $name]
		if {[regexp $nameSi $name]} {
			puts -nonewline $out [format "%-4s  " $typeSi]
			puts -nonewline $out [format "% 5.6f       " $chargeSi]
			puts -nonewline $out [format "%6.4f " $massSi]
		} else {
			puts -nonewline $out [format "%-4s  " $typeN]
			puts -nonewline $out [format "% 5.6f       " $chargeN]
			puts -nonewline $out [format "%6.4f " $massN]
		}
		puts $out $dummy
		
		incr atom
	}
	close $inStream
	puts $out ""
  
	##### BONDS
	# Write the bonds.
	set total [format %8i $totalBonds]
	puts $out "$total !NBOND: bonds"
	set num 0
	foreach b $bond {
		puts -nonewline $out [format "%8i%8i" [lindex $b 0] [lindex $b 1]]
	
		incr num
			if {$num == 4} {
			puts $out ""
			set num 0
		}
	}
	puts $out ""

	##### ANGLES
	# Write the angles.
	puts $out "[format %8i $totalAngles] !NTHETA: angles"
	set num 0
	foreach a $angle {
		puts -nonewline $out \
		[format "%8i%8i%8i" [lindex $a 0] [lindex $a 1] [lindex $a 2]]
	
		incr num
		if {$num == 3} {
			puts $out ""
			set num 0
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

	##### NON-BONDED
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



