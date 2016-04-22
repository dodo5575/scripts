# tclForce.tcl
# Inputs: targetAtomPdb forcesRecalcFreq dataFile
print "Starting tclForces: tclForce"

# Define the target beta value.
set targetMark "1.00"

set atomList {}
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

set outStream [open $dataFile w]
print "OUTFILE: $dataFile"
fconfigure $outStream -translation binary

###################################################################
# This procedure is executed at each time step.
set forcesRecalcFreq1 [expr $forcesRecalcFreq-1]
set recalcCount $forcesRecalcFreq1
print "Starting calcforces..."
proc calcforces {} {
    global recalcCount forcesRecalcFreq forcesRecalcFreq1
    global targetAtoms pos0 outStream

    if {$recalcCount == $forcesRecalcFreq1} {
	loadcoords pos0
    } elseif {$recalcCount >= $forcesRecalcFreq} {
	loadtotalforces force
	loadcoords pos
	
	foreach index $targetAtoms {
	    #print "TCLFORCES $step $index $force($index)"
	    #print "TCLFORCES $step $index"
	    set data [concat $pos0($index) $pos($index) $force($index)]
	    puts -nonewline $outStream [binary format f9 $data]
	    flush $outStream
	}
	set recalcCount 0
    }

    incr recalcCount
    return
}
