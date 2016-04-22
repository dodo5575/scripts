# Author: Jeff Comer <jcomer2@illinois.edu>
# Parameters:
set s0 2.5
set z0 25
set systemSize {40.0 40.0 72.0}
set dl 1.0
set range 8.0
# Input:
set inProfile surf_profile.dat
# Output:
set outGrid built_pmf_no.dx

source vector.tcl
source gridForce.tcl

proc readProfile {fileName} {
    set in [open $fileName r]
    set prof {}
    
    while {[gets $in line] >= 0} {
	if {[string length $line] <= 1} { continue }
	if {[string equal [string index $line 0] "\#"] } { continue }
	
	set tok [concat $line]
	set s [lindex $tok 0]
	set v [lindex $tok 1]
	lappend prof [list $s $v]
    }
    
    close $in
    return $prof
}

proc maxValue {items} {
    set maxX [lindex $items 0]
    foreach x $items {
	if {$x > $maxX} {
	    set maxX $x
	}
    }
    return $maxX
}

# Interpolation between p1 and p2.
proc interpolate {posX p0 p1 p2 p3} {
    for {set i 0} {$i < 4} {incr i} {
	set x($i) [lindex [set p${i}] 0]
	set f($i) [lindex [set p${i}] 1]
    }
    if {$x(2) - $x(1) == 0.0} { return $f(1) }
    set wx [expr ($posX - $x(1))/($x(2) - $x(1))]

    # Mix.
    set a3 [expr 0.5*(-$f(0) + 3*$f(1) - 3*$f(2) + $f(3))]
    set a2 [expr 0.5*(2*$f(0) - 5*$f(1) + 4*$f(2) - $f(3))]
    set a1 [expr 0.5*(-$f(0) + $f(2))]
    set a0 [expr $f(1)]

    return [expr $a3*$wx*$wx*$wx + $a2*$wx*$wx + $a1*$wx + $a0]  
}

set prof [readProfile $inProfile]
set profX {}
set profV {}
foreach pair $prof {
    lappend profX [lindex $pair 0]
    lappend profV [lindex $pair 1]
}
set profDx [expr [lindex $profX 1]- [lindex $profX 0]]
set profX0 [lindex $profX 0]
set profX1 [lindex $profX end]
set profX0a [expr $profX0 + 2.0*$profDx]
set profX1a [expr $profX1 - 2.0*$profDx]
set profN [llength $profX]
set maxV [maxValue $profV]
puts "Profile from x=${profX0} to x=${profX1} of $profN points with a maximum of $maxV."
set zs [expr $z0-$range]

foreach {lx ly lz} $systemSize { break }
newGridFit pot $lx $ly $lz $dl

# Loop through the grid points.
set v {}
for {set ix 0} {$ix < $pot(nx)} {incr ix} {
    for {set iy 0} {$iy < $pot(ny)} {incr iy} {
	for {set iz 0} {$iz < $pot(nz)} {incr iz} {
	    #set ind [expr $iz + $iy*$pot(nz) + $ix*$pos(nz)*$pot(ny)]
	    set r [getPosition pot [list $ix $iy $iz]]

	    foreach {x y z} $r { break }
	    set s [expr sqrt($x*$x + $y*$y)]

	    # Handle the the pore walls, the membrane surface, and in between.
	    if {abs($z) < $zs} {
		set x $s
	    } elseif {$s > $range} {
		set x [expr $z0-abs($z)]
	    } else {
		set dz [expr abs($z)-$zs]
		set ds [expr $range - $s]
		set x [expr $range - sqrt($dz*$dz + $ds*$ds)] 
	    }

	    if {$s < $s0 || abs($z) > $z0} {
		lappend v 0.0
	    } elseif {($s > $range && abs($z) < $zs)} {
		lappend v $maxV
	    } else {
		# Interpolate from the PMF profile.
		set ind [expr int(floor(($x-$profX0)/$profDx))]
       
		# Get the coordinates for interpolation.
		for {set i 0} {$i < 4} {incr i} {
		    set indI [expr $ind - 1 + $i]
	
		    if {$indI < 0} { set indI 0 }
		    if {$indI >= $profN} { set indI [expr $profN-1] }
		    set px [lindex $profX $indI]
		    set pv [lindex $profV $indI]

		    set p${i} [list $px $pv]
		}
		lappend v [interpolate $x $p0 $p1 $p2 $p3]
		#lappend v 1.0
		#lappend v [lindex $p1 1]
	    }
	}
    }
    puts "[expr (100*$ix)/$pot(nx)] percent complete."
}
# Write the new grid.
set pot(data) $v
writeDx pot $outGrid

