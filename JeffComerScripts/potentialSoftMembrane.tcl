# potentialScale.tcl
# Scale by a factor.
# Requires gridForce.tcl.
# Author: Jeff Comer <jcomer2@uiuc.edu>

if {$argc != 5} {
    puts "This script requires five arguments: "
    puts "  basis file"
    puts "  grid size"
    puts "  membrane length"
    puts "  softness length"
    puts "  output file"
    puts ""
    exit
}


source vector.tcl
source gridForce.tcl

# Input:
set basisFile [lindex $argv 0]
# Output:
set outFile [lindex $argv end]
# Parameters:
set dl [expr 1.0*[lindex $argv 1]]
set z0 [expr 0.5*[lindex $argv 2]]
set l [lindex $argv 3]

proc readBasis {fileName} {
    set in [open $fileName r]

    set basis {}
    foreach line [split [read $in] "\n"] {
	if {[string equal [string index $line 0] "\#"]} {continue}
	set tok [concat $line]
	
	if {[llength $tok] < 3} {continue}
	lappend basis [lrange $tok 0 2]
    }
    close $in
    
    return $basis
}

proc centerBasis {basis} {
    set rx [vecScale -0.5 [lindex $basis 0]]
    set ry [vecScale -0.5 [lindex $basis 1]]
    set rz [vecScale -0.5 [lindex $basis 2]]
    set ret [vecAdd $rx $ry]
    set ret [vecAdd $ret $rz]
    return $ret
}

proc subdivideBasis {basis nx ny nz} {
    set ex [vecScale [expr 1.0/$nx] [lindex $basis 0]]
    set ey [vecScale [expr 1.0/$ny] [lindex $basis 1]]
    set ez [vecScale [expr 1.0/$nz] [lindex $basis 2]]
    return [list $ex $ey $ez]
}

set basis [readBasis $basisFile]
foreach {bx by bz} $basis {break}
puts "Read the basis: $basis"
set nx [expr int(ceil([vecLength $bx]/$dl))]
set ny [expr int(ceil([vecLength $by]/$dl))]
set nz [expr int(ceil([vecLength $bz]/$dl))]
set origin [centerBasis $basis]
set delta [subdivideBasis $basis $nx $ny $nz]
puts "Subdivision: $nx $ny $nz"
puts "New basis: $delta"
newGrid grid $delta $nx $ny $nz

proc potential {r} {
    global z0 l

    foreach {x y z} $r {break}
    set p [expr abs($z)]

    if {$p < $z0} {
	return 0.0
    } elseif {$p > $z0 + 2.0*$l} {
	return 1.0
    } else {
	set q [expr ($p-($z0+$l))/$l]
	return [expr 0.5 + 0.75*$q - 0.25*$q*$q*$q]
    }
}

# Calculate the grid values.
set n $grid(size)
set pot {}
for {set j 0} {$j < $n} {incr j} {
    set r [getPosition grid [latticePosition grid $j]]
    lappend pot [expr -[potential $r]]
}
set grid(data) $pot

writeDx grid $outFile
puts "Wrote $outFile."


