# mapWham.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

#set name anneal
set binN 30
# Parameters:
set tol 0.0002
set waterN 1305
set springK 4.0; # in kcal mol^1 A^-2
set posXList {-12.5 -7.5 -2.5 2.5 7.5}
set posYList {-12.5 -7.5 -2.5 2.5 7.5}
set cutTime 200000
set farZ0 18.5
set farZ1 19.5
# Input:
set simFile phos_pos.txt
set inPrefix pmf_${name}/dmmp_${name}_pos
set inSuffix .log.force
# Output:
set outPrefix quest_${name}

if {[string equal $name anneal]} {
    set inPrefix pmf_${name}/dmmp_map_pos
}
source whamPmf.tcl

set kBT 0.5862292; # in kcal/mol (at 295 K)
set waterVol 29.922; # in A^3
set avogadroNum 6.02214179e23
set volume [expr $waterN*$waterVol]
set kappa [expr $springK/$kBT]; # spring constant in k_B*T/A^2

set posList {}
foreach x $posXList {
    foreach y $posYList {
	lappend posList [list $x $y]
    }
}
puts "Mapping [llength $posList] pmf values."

set outData [open $outPrefix.txt w]
set outDepth [open $outPrefix.depth w]

foreach pos $posList {
    foreach {posX posY} $pos {break}
    puts "\nposition: $posX $posY"
    puts $outData "\nposition: $posX $posY"

    # Read the index file.
    set m [readData $simFile]
    set simIndex {}
    set simPos {}
    set simZ {}
    foreach item $m {
	if {[lindex $item 1] == $posX && [lindex $item 2] == $posY} {
	    lappend simIndex [lindex $item 0]
	    lappend simPos [lrange $item 1 3]
	    lappend simZ [lindex $item 3]
	}
    }

    # Read the position data.
    set data {}
    foreach s $simIndex {
	set dataFile ${inPrefix}${s}${inSuffix}
	
	set r [readPositions $dataFile $cutTime]
	puts "Loaded [llength $r] data points from $dataFile."

	set z {}
	foreach p $r {
	    lappend z [lindex $p 2]
	}

	lappend data $z
    }

    # Run WHAM.
    set pmf [wham data $simZ $kappa $binN $tol]

    # Find the pmf far from the surface.
    set freeEnergy {}
    set binCen {}
    set farU 0.0
    set farN 0
    foreach item $pmf {
	foreach {c u} $item {break}
	lappend freeEnergy $u
	lappend binCen $c

	if {$c > $farZ0 && $c < $farZ1} {
	    set farU [expr $farU + $u]
	    incr farN
	}
    }
    set farU [expr $farU/$farN]

    # Write the pmf.
    set outFile ${outPrefix}_${posX}_${posY}.dat
    set out [open $outFile w]
    foreach item $pmf {
	foreach {c u} $item {break}

	puts $out "$c [expr $u-$farU]"
    }
    close $out

    set bindingEnergy $farU
    set bindingConst [expr 1/$volume*exp($bindingEnergy/$kBT)]; # in A^-3
    set bindingConstMolar [expr 1e18*$bindingConst/$avogadroNum]; # in nM
    set bindingAffinity [expr -log10($bindingConstMolar)]
    
    puts "binding energy: $farU k_B T"
    puts "binding constant: $bindingConstMolar nM"
    puts "binding affinity: $bindingAffinity"

    puts $outData "binding energy: $farU k_B T"
    puts $outData "binding constant: $bindingConstMolar nM"
    puts $outData "binding affinity: $bindingAffinity"

    puts $outDepth "$posX $posY $bindingEnergy"
}

close $outData
close $outDepth
