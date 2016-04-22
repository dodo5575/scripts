# gridSourcePoints.tcl
# Puts points on a grid given a pore geometry.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 6} {
    puts "This script requires six arguments: "
    puts " basis file"
    puts " grid size"
    puts " pore length"
    puts " pore diameter"
    puts " pore angle"
    puts " output file"
    puts ""
    exit
}

set dl [expr 1.0*[lindex $argv 1]]
set poreLength [lindex $argv 2]
set poreDiameter [lindex $argv 3]
set poreAngle [lindex $argv 4]
# Input:
set basisFile [lindex $argv 0]
# Output:
set outFile [lindex $argv end]

# Requires vector.tcl
source vector.tcl

set pi [expr 4.0*atan(1.0)]
set s0 [expr 0.5*$poreDiameter]
set z0 [expr 0.5*$poreLength]
set slope [expr tan($poreAngle*$pi/180.0)]

proc inSpace {r s0 z0 slope} {
    foreach {x y z} $r {break}
    
    if {abs($z) > $z0} {return 1}
    set s [expr $s0 + $slope*abs($z)]
    if {$x*$x + $y*$y > $s*$s} {return 0}
    return 1
}

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

# Write the source points.
puts "Sampling [expr $nx*$ny*$nz] source points..."
set n 0
set out [open $outFile w]
for {set ix 0} {$ix < $nx} {incr ix} {
    for {set iy 0} {$iy < $ny} {incr iy} {
	for {set iz 0} {$iz < $nz} {incr iz} {
	    set r [vecAdd [vecTransform $delta [list $ix $iy $iz]] $origin]

	    if {[inSpace $r $s0 $z0 $slope]} {
		puts $out $r
		incr n
	    }
	}
    }
}
close $out
puts "Wrote $n source points to $outFile."



