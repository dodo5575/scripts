#
# to use:  vmd -dispdev text  -e computeDNAcurrentCONSTRICTION-A3.tcl
#

# define trajectory and structure files

set psffile /home/alek/HEMOLYSIN/30_AMIT-58/MUTATE_58/EQUIL-C3/final-C3.psf
set pdbfile /home/alek/HEMOLYSIN/30_AMIT-58/MUTATE_58/EQUIL-C3/final-C3.psf

set dcdfile /Scr/alek/GRID_FORCE-COBALT/GRID_FORCE-V3-C/gridForce-C3_Runs0-1move2-7_25ps_wrap.dcd
#set dcdfile tmp.dcd

set frameRate  25000
set skipframes 1


# load trajectory and wait until loading completes 
set molTRJ [mol load psf $psffile dcd $dcdfile]
set n [ molinfo top get numframes ]
puts "  Loaded $n frames."

# for A nucleotides:
set nucleotideMass 312.203998448

# for C nucleotides:
# set nucleotideMass 288.178398968


set dcdfile ${dcdfile}.CONSTRICTION
set geometryFile channelFluctuationsCONSTRICTION.dat

############################################
set name CLA
set outfile ${dcdfile}.${name}
set outCh [open ${outfile} w]


set inCh1 [open $geometryFile r]

set i 0
set halfSpacing 1.5
set selection ""
foreach line [split [read $inCh1] \n] {
    if {[llength $line] == 3} {
	foreach {z R var} $line { break }

	if {$i == 0} {
	    set zDown [expr $z - $halfSpacing]
	    set selection "name $name and ((z > [expr $z -$halfSpacing] and z < [expr $z + $halfSpacing] and x*x +y*y < $R*$R)"
	} else {

	    set selection "$selection or (z > [expr $z -$halfSpacing] and z < [expr $z + $halfSpacing] and x*x +y*y < $R*$R)"
	}
	incr i	
    }
}
set selection "$selection )"
set zUp [expr $z + $halfSpacing]
puts "zUp: $zUp; zDonw: $zDown"
close $inCh1

puts $selection

# define first and last frames and stride
set firstframe 0
set lastFrame $n


######

# Main loop

if {$name == "CLA"} {
    set sign -1.0
} elseif {$name == "POT"} {
    set sign 1.0
} else {

    puts "name $name has not been recognized"
    puts "EXIT"
    exit
}

puts "zUp = $zUp; zDown = $zDown"
set scaleFactor [expr $sign*1.60276e-19*1.e+15/($zUp -$zDown)/($skipframes*$frameRate)]
puts "scaleFactor $scaleFactor"

set listCLA ""

for {set frame 0} {$frame < $lastFrame} {incr frame $skipframes} {
    
    
    if {$frame != 0} {

	set currSelection [atomselect top $selection frame $frame]

	set currCoord [$currSelection get {index z}]
	array unset coord
	set indexCurr {}
	foreach line $currCoord {
	    foreach {index z} $line { break }
	    lappend indexCurr $index
	    array set coord [list $index  $z]
	}

#	puts "frame = $frame [array size coord]  [array size refCoord]"

	set dn 0.0
	set nWaters 0
	foreach index $indexCurr {

	    if {[array get refCoord $index] != ""} {
#		puts "$coord($index) $refCoord($index)"
		set dn [expr $dn + ($coord($index) - $refCoord($index))]
		incr nWaters
	    } else {
		set tmpZ  $coord($index)
		set tmpSel [atomselect top "index $index" frame [expr $frame -1] ]
		set oldZ [$tmpSel get z]
		if {$oldZ > $zUp} { 
		    set dn [expr $dn + ($tmpZ -$zUp)] 
		} elseif {$oldZ < $zDown} {
		    set dn [expr $dn + ($tmpZ -$zDown)]
		} else {
		    set dn [expr $dn + ($tmpZ -$oldZ)]
		}
		incr nWaters
		$tmpSel delete		
	    }
	}

#	puts "Only in refCoord"
	foreach index $indexRef {
	    if {[array get coord $index] == ""} {
		set tmpZ $refCoord($index)
		set tmpSel [atomselect top "index $index" frame $frame]
		set newZ [$tmpSel get z]
		if {$newZ > $zUp} { 
		    set dn [expr $dn + ($zUp-$tmpZ)] 
		} elseif {$newZ < $zDown} {
		    set dn [expr $dn + ($zDown-$tmpZ)]
		} else {
		    set dn [expr $dn + ($newZ -$tmpZ)]
		}
		incr nWaters
		$tmpSel delete
	    }
	}


	set dn [expr $dn*$scaleFactor]

	puts "frame = $frame    dn = $dn $nWaters"
	puts $outCh "$frame $dn $nWaters"
	lappend listCLA "$frame $dn $nWaters"
	
	set indexRef {}
	array unset refCoord
	foreach line $currCoord {
	    foreach {index z} $line { break }
	    lappend indexRef $index
	    array set refCoord [list $index  $z]	    
	}

    } else {

	set currSelection [atomselect top $selection frame $frame]

	set currCoord [$currSelection get {index z}]
	set indexRef {}
	foreach line $currCoord {
	    foreach {index z} $line { break }
	    lappend indexRef $index
	    array set refCoord [list $index  $z]
	}
    }
    $currSelection delete
}
close $outCh





########### REPEAT FOR K ###############


set name POT
set outfile ${dcdfile}.${name}
set outCh [open ${outfile} w]


set inCh1 [open $geometryFile r]

set i 0
set halfSpacing 1.5
set selection ""
foreach line [split [read $inCh1] \n] {
    if {[llength $line] == 3} {
	foreach {z R var} $line { break }

	if {$i == 0} {
	    set zDown [expr $z - $halfSpacing]
	    set selection "name $name and ((z > [expr $z -$halfSpacing] and z < [expr $z + $halfSpacing] and x*x +y*y < $R*$R)"
	} else {

	    set selection "$selection or (z > [expr $z -$halfSpacing] and z < [expr $z + $halfSpacing] and x*x +y*y < $R*$R)"
	}
	incr i	
    }
}
set selection "$selection )"
set zUp [expr $z + $halfSpacing]
puts "zUp: $zUp; zDonw: $zDown"
close $inCh1

puts $selection

set firstframe 0
set lastFrame $n

######

if {$name == "CLA"} {
    set sign -1.0
} elseif {$name == "POT"} {
    set sign 1.0
} else {

    puts "name $name has not been recognized"
    puts "EXIT"
    exit
}

set scaleFactor [expr $sign*1.60276e-19*1.e+15/($zUp -$zDown)/($frameRate*$skipframes)]
puts "scaleFactor $scaleFactor"

set listPOT ""

for {set frame 0} {$frame < $lastFrame} {incr frame $skipframes} {
    
    
    if {$frame != 0} {

	set currSelection [atomselect top $selection frame $frame]

	set currCoord [$currSelection get {index z}]
	array unset coord
	set indexCurr {}
	foreach line $currCoord {
	    foreach {index z} $line { break }
	    lappend indexCurr $index
	    array set coord [list $index  $z]
	}

#	puts "frame = $frame [array size coord]  [array size refCoord]"

	set dn 0.0
	set nWaters 0
	foreach index $indexCurr {

	    if {[array get refCoord $index] != ""} {
#		puts "$coord($index) $refCoord($index)"
		set dn [expr $dn + ($coord($index) - $refCoord($index))]
		incr nWaters
	    } else {
		set tmpZ  $coord($index)
		set tmpSel [atomselect top "index $index" frame [expr $frame -1] ]
		set oldZ [$tmpSel get z]
		if {$oldZ > $zUp} { 
		    set dn [expr $dn + ($tmpZ -$zUp)] 
		} elseif {$oldZ < $zDown} {
		    set dn [expr $dn + ($tmpZ -$zDown)]
		} else {
		    set dn [expr $dn + ($tmpZ -$oldZ)]
		}
		incr nWaters
		$tmpSel delete		
	    }
	}

#	puts "Only in refCoord"
	foreach index $indexRef {
	    if {[array get coord $index] == ""} {
		set tmpZ $refCoord($index)
		set tmpSel [atomselect top "index $index" frame $frame]
		set newZ [$tmpSel get z]
		if {$newZ > $zUp} { 
		    set dn [expr $dn + ($zUp-$tmpZ)] 
		} elseif {$newZ < $zDown} {
		    set dn [expr $dn + ($zDown-$tmpZ)]
		} else {
		    set dn [expr $dn + ($newZ -$tmpZ)]
		}
		incr nWaters
		$tmpSel delete
	    }
	}


	set dn [expr $dn*$scaleFactor]

	puts "frame = $frame    dn = $dn  $nWaters"
	puts $outCh "$frame $dn $nWaters"
	
	lappend listPOT "$frame $dn $nWaters"

	set indexRef {}
	array unset refCoord
	foreach line $currCoord {
	    foreach {index z} $line { break }
	    lappend indexRef $index
	    array set refCoord [list $index  $z]	    
	}

    } else {

	set currSelection [atomselect top $selection frame $frame]

	set currCoord [$currSelection get {index z}]
	set indexRef {}
	foreach line $currCoord {
	    foreach {index z} $line { break }
	    lappend indexRef $index
	    array set refCoord [list $index  $z]
	}
    }
    $currSelection delete
}
close $outCh

###### END ########

set outfile ${dcdfile}.ALL
set outCh [open ${outfile} w]

foreach line1 $listCLA line2 $listPOT {

    foreach {frame dnCLA nCLA} $line1 { break }
    foreach {frame dnPOT nPOT} $line2 { break }

    puts $outCh "$frame [expr $dnCLA + $dnPOT] [expr $nCLA + $nPOT]"
    
}
close $outCh




########### REPEAT FOR DNA ###############


set name DNA
set outfile ${dcdfile}.${name}
set outCh [open ${outfile} w]



set inCh1 [open $geometryFile r]

set i 0
set halfSpacing 1.5
set selection ""
foreach line [split [read $inCh1] \n] {
    if {[llength $line] == 3} {
	foreach {z R var} $line { break }

	if {$i == 0} {
	    set zDown [expr $z - $halfSpacing]
	    set selection "$name and ((z > [expr $z -$halfSpacing] and z < [expr $z + $halfSpacing] and x*x +y*y < $R*$R)"
	} else {

	    set selection "$selection or (z > [expr $z -$halfSpacing] and z < [expr $z + $halfSpacing] and x*x +y*y < $R*$R)"
	}
	incr i	
    }
}
set selection "$selection )"
set zUp [expr $z + $halfSpacing]
puts "zUp: $zUp; zDonw: $zDown"
close $inCh1

puts $selection

set firstframe 0
set lastFrame $n

######

if {$name == "CLA"} {
    set sign -1.0
} elseif {$name == "POT" || $name == "DNA" } {
    set sign 1.0
} else {

    puts "name $name has not been recognized"
    puts "EXIT"
    exit
}

# Units are: nucleotides / ns

set scaleFactor [expr $sign*1.e+6/($zUp -$zDown)/($frameRate*$skipframes)/$nucleotideMass ]
puts "scaleFactor $scaleFactor"

set listDNA ""

for {set frame 0} {$frame < $lastFrame} {incr frame $skipframes} {
    
    
    if {$frame != 0} {

	set currSelection [atomselect top $selection frame $frame]

	set currCoord [$currSelection get {index z mass}]
	array unset coord
	set indexCurr {}
	set indexMass {}
	foreach line $currCoord {
	    foreach {index z mass} $line { break }
	    lappend indexCurr $index
	    lappend indexMass $mass
	    array set coord [list $index  $z]
	}

#	puts "frame = $frame [array size coord]  [array size refCoord]"

	set dn 0.0
	set nWaters 0
	foreach index $indexCurr mass $indexMass {

	    if {[array get refCoord $index] != ""} {
#		puts "$coord($index) $refCoord($index)"
		set dn [expr $dn + $mass*($coord($index) - $refCoord($index))]
		incr nWaters
	    } else {
		set tmpZ  $coord($index)
		set tmpSel [atomselect top "index $index" frame [expr $frame -1] ]
		set oldZ [$tmpSel get z]
		if {$oldZ > $zUp} { 
		    set dn [expr $dn + $mass*($tmpZ -$zUp)] 
		} elseif {$oldZ < $zDown} {
		    set dn [expr $dn + $mass*($tmpZ -$zDown)]
		} else {
		    set dn [expr $dn + $mass*($tmpZ -$oldZ)]
		}
		incr nWaters
		$tmpSel delete		
	    }
	}

#	puts "Only in refCoord"
	foreach index $indexRef mass $massRef {
	    if {[array get coord $index] == ""} {
		set tmpZ $refCoord($index)
		set tmpSel [atomselect top "index $index" frame $frame]
		set newZ [$tmpSel get z]
		if {$newZ > $zUp} { 
		    set dn [expr $dn + $mass*($zUp-$tmpZ)] 
		} elseif {$newZ < $zDown} {
		    set dn [expr $dn + $mass*($zDown-$tmpZ)]
		} else {
		    set dn [expr $dn + $mass*($newZ -$tmpZ)]
		}
		incr nWaters
		$tmpSel delete
	    }
	}


	set dn [expr $dn*$scaleFactor]

	puts "frame = $frame    dn = $dn  $nWaters"
	puts $outCh "$frame $dn $nWaters"
	
	lappend listDNA "$frame $dn $nWaters"

	set indexRef {}
	set massRef {}
	array unset refCoord
	foreach line $currCoord {
	    foreach {index z mass} $line { break }
	    lappend indexRef $index
	    lappend massRef $mass
	    array set refCoord [list $index  $z]	    
	}

    } else {

	set currSelection [atomselect top $selection frame $frame]

	set currCoord [$currSelection get {index z mass}]
	set indexRef {}
	set massRef {}
	foreach line $currCoord {
	    foreach {index z mass} $line { break }
	    lappend indexRef $index
	    lappend massRef $mass
	    array set refCoord [list $index  $z]
	}
    }
    $currSelection delete
}
close $outCh

###### END ########




quit
