set poreLength 50
set poreDiameter 11
set poreAngle 10
set dl 0.5
set dx 40.0
set dy 40.0
set dz 60.0
# Output:
set outFile phantom0.5.txt

set pi [expr 4.0*atan(1.0)]
set s0 [expr 0.5*$poreDiameter]
set z0 [expr 0.5*$poreLength]
set slope [expr tan($poreAngle*$pi/180.0)]

set nx [expr int(ceil($dx/$dl))]
set ny [expr int(ceil($dy/$dl))]
set nz [expr int(ceil($dz/$dl))]

set xi [expr -0.5*$dx]
set yi [expr -0.5*$dy]
set zi [expr -0.5*$dz]

proc inPore {x y z} {
    global s0 z0 slope
    
    set rad [expr $s0 + $slope*abs($z)]
    return [expr $x*$x + $y*$y > $rad*$rad && abs($z) < $z0]
}

puts "Contains [expr $nx*$ny*$nz] points!"

set out [open $outFile w]
for {set ix 0} {$ix < $nx} {incr ix} {
    for {set iy 0} {$iy < $ny} {incr iy} {
	for {set iz 0} {$iz < $nz} {incr iz} {
	    set x [expr $xi + $dx*$ix/$nx]
	    set y [expr $yi + $dy*$iy/$ny]
	    set z [expr $zi + $dz*$iz/$nz]
	    
	    if {[inPore $x $y $z]} {
		puts $out "$x $y $z"
	    }
	}
    }
}

close $out



