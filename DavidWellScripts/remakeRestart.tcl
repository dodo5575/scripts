source $env(SCRIPTS)/Procs.tcl

package require math::statistics

vmdargs remakeRestart.tcl psf dcd dcdrate temperature outbase

# figure out thermal velocity distribution
proc thermal { mol } {
    global temperature
    set kT [expr {$temperature * 0.0019858775}]
    set all [atomselect $mol all]
    set n [$all num]
    $all delete
    set vels {}
    for { set ind 0 } { $ind < $n } { incr ind } {
	progressbar $ind $n 40
	set sel [atomselect $mol "index $ind"]
	set mass [lindex [$sel get mass] 0]
	set stddev [expr {sqrt($kT / $mass)}]
	set vel [::math::statistics::random-normal 0 $stddev 3]
	push vels $vel
	$sel delete
    }
    return $vels
}

regsub "dcd$" $dcd "xst" xst
regsub "dcd$" $dcd "info.txt" info

set frame [expr [exec totalframes $dcd] - 1]

mol load psf $psf
mol addfile $dcd first $frame last $frame waitfor all

# write coor
set all [atomselect top all]
$all writenamdbin $outbase.coor
puts "\nWROTE $outbase.coor\n"

# write vel
$all set {x y z} [thermal top]
$all writenamdbin $outbase.vel
puts "\nWROTE $outbase.vel\n"

# write xsc
set inch [open $xst r]
set flag 0
foreach line [split [read $inch] "\n"] {
    if { [regexp "^#" $line] } {
	continue
    }
    
    lassign [split $line] ts ax ay az bx by bz cx cy cz ox oy oz
    
    if { ! $flag } {
	set flag 1
	set step [expr $dcdrate * $frame - $ts]
	puts "step = $step"
    }
    if { $ts == $step } {
	set ouch [open $outbase.xsc w]
	puts $ouch "# NAMD extended system configuration restart file"
	puts $ouch "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z"
	puts $ouch $line
	
	set dims_dcd [molinfo top get {a b c}]
	set dims_xst [list $ax $by $cz]
	puts "delta: [vecsub $dims_dcd $dims_xst]"
	close $ouch
	puts "\nWROTE $outbase.xsc\n"
	break
    }
}

close $inch

set ouch [open $info w]
puts $ouch "# created with $env(SCRIPTS)/remakeRestart.tcl"
fputvar $ouch psf dcd dcdrate temperature outbase
close $ouch

exit
