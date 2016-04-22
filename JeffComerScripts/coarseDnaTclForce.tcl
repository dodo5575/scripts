######################################################################
# Author: Jeff Comer <jcomer2@illinois.edu>
######################################################################
# Inputs: coarseDnaPdb, stackingPairsFile, ionConcentration, elecTemperature
print "Starting coarseDnaTclForce"
print "coarseDnaPdb: $coarseDnaPdb"
print "stackingPairsFile: $stackingPairsFile"

set temp 295
set ionicConc 0.150; # in mol/kg

set eps 0.26; # in kcal/mol
set epsCytGua [expr 4.0*$eps]
set epsAdeThy [expr 2.0/3.0*$epsCytGua]

set qPho -1.0; # phosphate charge (e)
set dielectricConst 2.398e-4*78.0; # dielectric const in water (mol e^2)/(kcal A)
set kCoulomb 1.0/(16.0*atan(1.0)*$dielectricConst); # Coulomb constant (kcal A/(mol e^2))
set kDebye 13.603; # Debye length (A)

set dCut 6.86
set sigMismatch [expr pow(2.0,-1.6)*1.0]
set sig0 [expr pow(2.0,-1.6)*$dCut]
set sigAdeThy 2.9002
set sigCytGua 2.8694

set targetMark 1.00

# Precompute powers.
set qPho2 [expr $qPho*$qPho]
set sigAdeThy12 [expr pow($sigAdeThy,12.0)]
set sigAdeThy10 [expr pow($sigAdeThy,10.0)]
set sigCytGua12 [expr pow($sigCytGua,12.0)]
set sigCytGua10 [expr pow($sigCytGua,10.0)]
set sig06 [expr pow($sig0,6.0)]
set sigMismatch6 [expr pow($sigMismatch,6.0)]

set outCharge [open pairlist_charge.txt w]
set outBase [open pairlist_base.txt w]
set outMismatch [open pairlist_mismatch.txt w]
set outEx [open pairlist_ex.txt w]

proc readNativeContacts {inFile} {
    set stackList {}
    set inStream [open $inFile r]
    foreach line [split [read $inStream] \n] {
	set tok [concat $line]
	if {[llength $tok] < 7} {continue}
	foreach {segName0 resId0 name0 segName1 resId1 name1 sig} $tok {break}
    
	set i0 [atomid $segName0 $resId0 $name0]
	set i1 [atomid $segName1 $resId1 $name1]
	set sig6 [expr pow($sig,6.0)]
	# Swap so that the smallest index is first.
	if {$i1 < $i0} {
	    set ind $i0
	    set i0 $i1
	    set i1 $ind
	}

	lappend stackList [list $i0 $i1 $sig6]
    } 
    close $inStream
    return $stackList
}

proc readAtoms {inFile} {
    global targetMark
    set atomList {}

    set inStream [open $inFile r]
    foreach line [split [read $inStream] \n] {
	if {![string match "ATOM*" $line]} {continue}
	
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

	if {$beta == $targetMark} {
	    lappend atomList [list $segName $resId $name]
	}
    }
    return $atomList
}

proc isChargePair {id0 id1} {
    foreach {s0 r0 n0} $id0 {s1 r1 n1} $id1 {break}
    
    if {![string equal $n0 PHO] || ![string equal $n1 PHO]} {
	return 0
    } elseif {[string equal $s0 $s1] && [expr abs($r1-$r0) <= 1]} {
	return 0
    } else {
	return 1
    }
}

proc isBasePair {id0 id1} {
    foreach {s0 r0 n0} $id0 {s1 r1 n1} $id1 {break}
    if {[string equal AB $n0] && [string equal TB $n1]} {return 1}
    if {[string equal AB $n1] && [string equal TB $n0]} {return 1}
    if {[string equal CB $n0] && [string equal GB $n1]} {return 2}
    if {[string equal CB $n1] && [string equal GB $n0]} {return 2}
    return 0
}

proc isMismatchPair {id0 id1} {
    foreach {s0 r0 n0} $id0 {s1 r1 n1} $id1 {break}
    if {[string equal AB $n0] && [string equal AB $n1]} {return 1}
    if {[string equal CB $n0] && [string equal CB $n1]} {return 1}
    if {[string equal GB $n0] && [string equal GB $n1]} {return 1}
    if {[string equal TB $n0] && [string equal TB $n1]} {return 1}

    if {[string equal AB $n0] && [string equal CB $n1]} {return 1}
    if {[string equal AB $n1] && [string equal CB $n0]} {return 1}
    if {[string equal AB $n0] && [string equal GB $n1]} {return 1}
    if {[string equal AB $n1] && [string equal GB $n0]} {return 1}
    if {[string equal CB $n0] && [string equal TB $n1]} {return 1}
    if {[string equal CB $n1] && [string equal TB $n0]} {return 1}
    if {[string equal GB $n0] && [string equal TB $n1]} {return 1}
    if {[string equal GB $n1] && [string equal TB $n0]} {return 1}
    return 0
}

proc isNonbondedPair {id0 id1} {
    foreach {s0 r0 n0} $id0 {s1 r1 n1} $id1 {break}
    
    # Check if they are within the same residue.
    if {$r0 == $r1} {
	# Check for bonds.
	if {[regexp ".B" $n0] && [string equal SUG $n1]} {return 0}
	if {[regexp ".B" $n1] && [string equal SUG $n0]} {return 0}
	if {[string equal PHO $n0] && [string equal SUG $n1]} {return 0}
	if {[string equal PHO $n1] && [string equal SUG $n0]} {return 0}
    }
    
    # Check if they are in adjacent residues.
    if {$r0+1 == $r1} {
	# Check for bonds.
	if {[string equal SUG $n0] && [string equal PHO $n1]} {return 0}
    }
    if {$r0-1 == $r1} {
	# Check for bonds.
	if {[string equal PHO $n0] && [string equal SUG $n1]} {return 0}
    }
    return 1
}

print "Reading DNA beads from $coarseDnaPdb..."
set beadList [readAtoms $coarseDnaPdb]
set nBeads [llength $beadList]
print "Registered $nBeads DNA beads."
print "There are [expr $nBeads*($nBeads-1)/2] pairs."

print "Extracting the bead indices."
set indexList {}
foreach bead $beadList {
    foreach {segName resId name} $bead {break}
    lappend indexList [atomid $segName $resId $name]
}

print "Reading stacking interactions from $stackingPairsFile..."
set stackList1 [readNativeContacts $stackingPairsFile]
set stackList {}
set stackSig6List {}
foreach stack $stackList1 {
    lappend stackList [list [lindex $stack 0] [lindex $stack 1]]
    lappend stackSig6List [lindex $stack 2]
}
print "Registered [llength $stackList] stacking pairs."

print "Generating pair lists..."
set chargeList {}
set baseList {}
set baseTypeList {}
set mismatchList {}
set exList {}
set nBeads1 [expr $nBeads-1]
for {set i 0} {$i < $nBeads1} {incr i} {
    for {set j [expr $i+1]} {$j < $nBeads} {incr j} {
	set idI [lindex $beadList $i]
	set idJ [lindex $beadList $j]
	set indexI [lindex $indexList $i]
	set indexJ [lindex $indexList $j]

	# Swap so that the smallest index is first.
	if {$indexJ < $indexI} {
	    set ind $indexI
	    set indexI $indexJ
	    set indexJ $ind
	}

	# Add valid charges to the charge pair list.
	if {[isChargePair $idI $idJ]} {
	    lappend chargeList [list $indexI $indexJ]
	    puts $outCharge "$idI $idJ"
	}

	# Don't add to the base, mismatch, or excluded volume pair list
	# if already a stacking pair.
	if {[lsearch $stackList [list $indexI $indexJ]] != -1} {continue}

	# Check nonbonded pairs for either base pair or excluded volume
	# interactions.
	if {[isNonbondedPair $idI $idJ]} {
	    set bpType [isBasePair $idI $idJ]
	    if {$bpType > 0} {
		lappend baseList [list $indexI $indexJ]
		lappend baseTypeList $bpType
		puts $outBase "$idI $idJ"
	    } elseif {[isMismatchPair $idI $idJ]} {
		lappend mismatchList [list $indexI $indexJ]
		puts $outMismatch "$idI $idJ"
	    } else {
		lappend exList [list $indexI $indexJ]
		puts $outEx "$idI $idJ"
	    }
	}
    }
}
close $outCharge
close $outBase
close $outMismatch
close $outEx

print "Registered [llength $chargeList] charge pairs."
print "Registered [llength $baseList] base pairs."
print "Registered [llength $mismatchList] mismatch base pairs."
print "Registered [llength $exList] excluded volume pairs."

# Add all of the DNA beads.
foreach i $indexList {
    addatom $i
}

###################################################################
# This procedure is executed at each time step.
print "Starting calcforces..."
proc calcforces {} {
    global stackList stackSig6List chargeList baseList baseTypeList mismatchList exList
    global qPho2 kCoulomb kDebye
    global eps epsAdeThy epsCytGua sigAdeThy10 sigAdeThy12 sigCytGua10 sigCytGua12
    global sigMismatch6 sig06 dCut

    loadcoords pos

    # Compute the charge interactions.
    foreach pair $chargeList {
	foreach {i j} $pair {break}
	set ri $pos($i)
	set rj $pos($j)

	set d [vecsub $ri $rj]
	set d1 [getbond $ri $rj]
	set d2 [expr $d1*$d1]

	set fScale [expr $qPho2*$kCoulomb/$d2*(1.0/$d1 + 1.0/$kDebye)*exp(-$d1/$kDebye)]
	set fi [vecscale $fScale $d]
	set fj [vecscale -1.0 $fi]
	addforce $i $fi
	addforce $j $fj

	addenergy [expr $qPho2*$kCoulomb/$d1*exp(-$d1/$kDebye)]
    }

    # Compute the stacking interactions.
    foreach pair $stackList sig6 $stackSig6List {
	foreach {i j} $pair {break}
	set ri $pos($i)
	set rj $pos($j)

	set d [vecsub $ri $rj]
	set d1 [getbond $ri $rj]
	set d8 [expr pow($d1,8)]
	set d14 [expr pow($d1,14)]

	set fScale [expr 4.0*$eps*(12.0*$sig6*$sig6/$d14 - 6.0*$sig6/$d8)]
	set fi [vecscale $fScale $d]
	set fj [vecscale -1.0 $fi]
	addforce $i $fi
	addforce $j $fj

	addenergy [expr 4.0*$eps*($sig6*$sig6/$d14 - $sig6/$d8)*$d1*$d1]
    }

    # Compute the base interaction.
    foreach pair $baseList type $baseTypeList {
	foreach {i j} $pair {break}
	set ri $pos($i)
	set rj $pos($j)

	set d [vecsub $ri $rj]
	set d1 [getbond $ri $rj]
	set d12 [expr pow($d1,12)]
	set d14 [expr pow($d1,14)]

	if {$type == 1} {
	    set epsBase $epsAdeThy
	    set sigBase10 $sigAdeThy10
	    set sigBase12 $sigAdeThy12
	} else {
	    set epsBase $epsCytGua
	    set sigBase10 $sigCytGua10
	    set sigBase12 $sigCytGua12
	}

	set fScale [expr 240.0*$epsBase*($sigBase12/$d14 - $sigBase10/$d12)]
	set fi [vecscale $fScale $d]
	set fj [vecscale -1.0 $fi]
	addforce $i $fi
	addforce $j $fj

	addenergy [expr 4.0*$epsBase*(5.0*$sigBase12/$d14 - 6.0*$sigBase10/$d12)*$d1*$d1]
    }
    
    # Compute the mismatch excluded volume interaction.
    foreach pair $mismatchList {
	foreach {i j} $pair {break}
	set ri $pos($i)
	set rj $pos($j)
	
	set d1 [getbond $ri $rj]
	if {$d1 > $dCut} {continue}
	
	set d8 [expr pow($d1,8)]
	set d14 [expr pow($d1,14)]
	set d [vecsub $ri $rj]

	set fScale [expr 4.0*$eps*(12.0*$sigMismatch6*$sigMismatch6/$d14 - 6.0*$sigMismatch6/$d8)]
	set fi [vecscale $fScale $d]
	set fj [vecscale -1.0 $fi]
	addforce $i $fi
	addforce $j $fj

	addenergy [expr 4.0*$eps*($sigMismatch6*$sigMismatch6/$d14 - $sigMismatch6/$d8)*$d1*$d1 + $eps]
    }

    # Compute the excluded volume interaction.
    foreach pair $exList {
	foreach {i j} $pair {break}
	set ri $pos($i)
	set rj $pos($j)
	
	set d1 [getbond $ri $rj]
	if {$d1 > $dCut} {continue}
	
	set d8 [expr pow($d1,8)]
	set d14 [expr pow($d1,14)]
	set d [vecsub $ri $rj]
	
	set fScale [expr 4.0*$eps*(12.0*$sig06*$sig06/$d14 - 6.0*$sig06/$d8)]
	set fi [vecscale $fScale $d]
	set fj [vecscale -1.0 $fi]
	addforce $i $fi
	addforce $j $fj

	addenergy [expr 4.0*$eps*($sig06*$sig06/$d14 - $sig06/$d8)*$d1*$d1 + $eps]
    }
}



