# softPore.tcl
# Force a selection out of double cone pore in a membrane.
# Inputs: springPore lengthPore radius0 radius1 forcesRecalcFreq targetAtomPdb
print "Starting tclForces: softPore"
set membraneMargin 4.0

set pi [expr 4.*atan(1.)]
set ds [expr $radius1-$radius0]
set slope [expr 2.0*$ds/$lengthPore]
set normalMag [expr sqrt($ds*$ds + 0.25*$lengthPore*$lengthPore)]
set normalS [expr -0.5*$lengthPore/$normalMag]
set normalZ [expr $ds/$normalMag]
# Define the target beta value.
set targetMark "1.00"

set atomList {}
set chargeList {}
print "Reading the list of target atoms..."

# The beta column determines whether the force should be exerted.
set inStream [open $targetAtomPdb r]
foreach line [split [read $inStream] \n] {
	# Extract each pdb field.
    set record [string trim [string range $line 0 5]]
		set serial [string trim [string range $line 6 10]]
		set name [string trim [string range $line 12 15]]
		set altLoc [string trim [string range $line 16 16]]
		set resName [string trim [string range $line 17 19]]
		set chainId [string trim [string range $line 21 21]]
		set resId [string trim [string range $line 22 25]]
		set iCode [string trim [string range $line 26 26]]
		set x [string trim [string range $line 30 37]]
		set y [string trim [string range $line 38 45]]
		set z [string trim [string range $line 46 53]]
		set occupancy [string trim [string range $line 54 59]]
		set beta [string trim [string range $line 60 65]]
		set segName [string trim [string range $line 72 75]]
		set element [string trim [string range $line 76 77]]
		set charge [string trim [string range $line 78 79]]

    if {[string match "ATOM*" $record] || [string match "HETA*" $record]} {
    		if {[string trim $beta] == $targetMark} {
			lappend atomList [list $segName $resId $name]
		}
    }  

}
close $inStream

# Make a list of indices.
set targetAtoms {}
foreach atom $atomList {
    foreach {segName resId name} $atom {break}
    set index [atomid $segName $resId $name]
    addatom $index
    lappend targetAtoms $index 
    print "TCLFORCESATOM $index $segName $resId $name"
}

if {[llength $targetAtoms] == 0} {
	print "WARNING: No atoms to which to apply tclForces."
}
print "Registered [llength $targetAtoms] target atoms."

# Return the force if the atom is within the pore, a scalar zero otherwise. 
proc poreForce {x y z} {   
    global springPore radius0 slope normalS normalZ lengthPore membraneMargin
    
    # Quit if we're above or below the pore.
    if {[expr $z < 0.0]} {return 0}
    if {[expr $z > 0.5*$lengthPore+$membraneMargin]} {return 0}

    # Compute the pore radius at this height.
    set s [expr sqrt($x*$x + $y*$y)]
    set sPore [expr $radius0 + $slope*$z]
    
    # Quit if we're inside the pore.
    if {[expr $s < $sPore]} {return 0}

    # Check for membrane or pore collision.
    if {[expr abs($z) > 0.5*$lengthPore]} {
	# Membrane collision.
	set normal [list 0 0 1]
	set d [expr 0.5*$lengthPore+$membraneMargin-$z]
    } else {
	# Pore collision.
	set normal [list [expr $normalS*$x/$s] [expr $normalS*$y/$s] [expr $normalZ]]
	set d [expr ($s-$radius0)*$normalS + $z*$normalZ]
    }

    # Apply a harmonic force away from the pore.
    return [vecscale [expr -$springPore*$d] $normal]
}

# Return the depth into the walls of the pore/membrane.
# Intended for diagnostics, such as determining springPore.
proc depth {x y z} {   
    global radius0 slope normalS normalZ lengthPore membraneMargin

    # Compute the pore radius at this height.
    set s [expr sqrt($x*$x + $y*$y)]
    
    # Check for membrane or pore collision.
    if {[expr abs($z) > 0.5*$lengthPore]} {
	# Membrane collision.
	if {[expr $z > 0.0]} {
	    set d [expr 0.5*$lengthPore+$membraneMargin-$z]
	} else {
	    set d [expr $z-0.5*$lengthPore-$membraneMargin]
	}
    } else {
	# Pore collision.
	set d [expr ($s-$radius0)*$normalS + $z*$normalZ]
    }

    return $d
}

###################################################################
# This procedure is executed at each time step.
set forceAtoms {}
set forces {}
set forcesRecalcFreq1 [expr $forcesRecalcFreq-1]
set recalcCount forcesRecalcFreq
print "Starting calcforces..."
print "Force sums: step number_of_atoms force_sum"
print "Forces: step index x y z Fx Fy Fz"
proc calcforces {} {
	global recalcCount forcesRecalcFreq forcesRecalcFreq1
	global targetAtoms pi
	global forceAtoms forces
		
	if {$recalcCount == $forcesRecalcFreq1} {
		# Get the atom coordinates for the next force calculation.
		foreach index $targetAtoms {
			addatom $index
		}
	} elseif {$recalcCount >= $forcesRecalcFreq} {
		# Recalculate the forces.
		set forceSum [list 0.0 0.0 0.0]
		loadcoords coords
		set step [getstep]
				
		set forceAtoms {}
		set forces {}
		foreach index $targetAtoms {
			foreach {x y z} $coords($index) {break}
			# Compute the force, if there is one.
			set f [poreForce $x $y $z]
			
			if {[llength $f] != 1} {
				lappend forceAtoms $index
				lappend forces $f
				
				print "TCLFORCE $step $index $x $y $z $f"
			        print "TCLFORCEDEPTH $step $index [depth $x $y $z]"
			        set forceSum [vecadd $forceSum $f]
			}
		}
		
		print "TCLFORCESSUM $step [llength $forces] $forceSum"
		clearconfig
		set recalcCount 0
	}
	
	# Add the forces.
	foreach i $forceAtoms f $forces {
	    addforce $i $f
	}
	
	incr recalcCount
	return
}



