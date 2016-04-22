######################################################################
# this is a tcl script for applying an oscillating electric field to
# selected atoms
# use with name tclforces
# written by alek@ks.uiuc.edu
# Author: Jeff Comer <jcomer2@illinois.edu>
######################################################################
# Inputs: switchZ, forcesRecalcFreq, targetAtomPdb

print "Starting tclElectricFieldSwitch"

# define all protein atoms :
set targetMark    "1.00"

set markList {}
set atomList {}
set chargeList {}
print "Reading the list of target atoms..."

# The beta column determines whether the force should be exerted and the
# occupancy column must contain the charge of the atom.
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

    if { [string match "ATOM*" $record] } {
	if {[string trim $beta] == $targetMark} {
    	    print "TCLFORCESMARK $segName $resId $name"
	    lappend markList [atomid $segName $resId $name]
	}
    }
}
close $inStream

# Check counts of atoms.
if {[llength $markList] == 0} {
    print "WARNING: No marked atoms for electric field switching!"
}
print "Registered [llength $markList] marked atoms."

###################################################################
# This procedure is executed at each time step.
set forcesRecalcFreq1 [expr $forcesRecalcFreq-1]
set recalcCount $forcesRecalcFreq1
set fieldOn 1
print "Starting calcforces..."
proc calcforces {} {
    global recalcCount forcesRecalcFreq forcesRecalcFreq1 switchZ
    global atomList forceList markList fieldOn markGroup

    if {$fieldOn} {
	if {$recalcCount == $forcesRecalcFreq1} {
	    # Get the atom coordinates for the next force calculation.
	    #set markGroup [addgroup $markList]
	    addatom [lindex $markList 0]
	    addatom [lindex $markList 1]
	} elseif {$recalcCount == $forcesRecalcFreq} {
	    # Get coordinates.
	    loadcoords coor
	    #set z [lindex $coor($markGroup) 2]
	    set z0 [lindex $coor([lindex $markList 0]) 2]
	    set z1 [lindex $coor([lindex $markList 1]) 2]
	    set z [expr {0.5*($z0+$z1)}]
	    print "TCLFORCESPOS $z"

	    # If the center of mass of the marked group is below switchZ, turn off the field.
	    if {$z < $switchZ} {
		set fieldOn 0
		print "Switching tclforces off!"
		print "TCLFORCESOFF [getstep]"
		eField 0 0 0
		clearconfig
		return
	    }

	    # Vital to use clearconfig after calculating forces, otherwise
	    # coordinates will be retrieved every timestep
	    clearconfig
	    set recalcCount 0
	}
    }
    
    incr recalcCount
    return    
}
