# Generate a periodic box from an xsc file.
# Author: Jeff Comer <jcomer2@illinois.edu>

proc periodicCell {xscFile {molid top}} {
    set pi [expr 4.0*atan(1.0)]
    set in [open $xscFile r]
    
    foreach line [split [read $in] "\n"] {
	if {[llength $line] < 10} {continue}
	if {[string match "#*" $line]} {continue}
	
	set tok [concat $line]
	foreach {step ax ay az bx by bz cx cy cz ox oy oz} $tok {break}
	
	set aMag [expr sqrt($ax*$ax + $ay*$ay + $az*$az)]
	set bMag [expr sqrt($bx*$bx + $by*$by + $bz*$bz)]
	set cMag [expr sqrt($cx*$cx + $cy*$cy + $cz*$cz)]
	
	set b_c [expr $bx*$cx + $by*$cy + $bz*$cz]
	set a_c [expr $ax*$cx + $ay*$cy + $az*$cz]
	set a_b [expr $ax*$bx + $ay*$by + $az*$bz]
	
	set alpha [expr acos($b_c/$bMag/$cMag)*180./$pi]
	set beta [expr acos($a_c/$aMag/$cMag)*180./$pi]
	set gamma [expr acos($a_b/$aMag/$bMag)*180./$pi]
    }
    close $in
    
    puts "alpha $alpha"
    puts "beta $beta"
    puts "gamma $gamma"
    puts "a $aMag"
    puts "b $bMag"
    puts "c $cMag"
    
    if {[string equal $molid all]} {
	set molList [molinfo list]
    } else {
	set molList [list $molid]
    }

    foreach m $molList {
	molinfo $m set alpha $alpha
	molinfo $m set beta $beta
	molinfo $m set gamma $gamma
	molinfo $m set a $aMag
	molinfo $m set b $bMag
	molinfo $m set c $cMag
    }

    return [list $aMag $bMag $cMag $alpha $beta $gamma]
}

proc periodicCellMan {lx ly lz {molid top}} {
    set pi [expr 4.0*atan(1.0)]
    puts "a $lx"
    puts "b $ly"
    puts "c $lz"
    
    if {[string equal $molid all]} {
	set molList [molinfo list]
    } else {
	set molList [list $molid]
    }

    foreach m $molList {
	molinfo $m set alpha 90
	molinfo $m set beta 90
	molinfo $m set gamma 90
	molinfo $m set a $lx
	molinfo $m set b $ly
	molinfo $m set c $lz
    }
}

proc setCell {a b c alpha beta gamma {molid top}} {
    molinfo $molid set alpha $alpha
    molinfo $molid set beta $beta
    molinfo $molid set gamma $gamma
    molinfo $molid set a $a
    molinfo $molid set b $b
    molinfo $molid set c $c
}

proc setCellFrames {a b c alpha beta gamma {molid top}} {
    set numFrames [molinfo $molid get numframes]
    for {set i 0} {$i < $numFrames} {incr i} {
	molinfo $molid set frame $i
	setCell $a $b $c $alpha $beta $gamma
    }
}