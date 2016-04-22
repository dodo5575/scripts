# Author: jcomer2@illinois.edu

set segA ADNA
set resA 4
set segB BDNA
set resB 108
set nBasepairs 4
# Input:
set pmfGrid ../wham/pmf_pmf_at_helix_pot.dx
set localGrid basepair_zero.dx
set psf ../one_turn_at_pot.psf
set pdb ../mark_dna_at_pot.pdb
# Output:
set outGrid pmf_local_at_helix1.dx

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

proc getBasepairTransform {segA resA segB resB {dnaMol top}} {
    set basis [getBasepairBasisBest $segA $resA $segB $resB $dnaMol]
    set origin [getBasepairPos $segA $resA $segB $resB $dnaMol]
    return [list $basis $origin]
}

proc extractBasepairPmf {} {
    global pmfGrid localGrid segA resA segB resB nBasepairs psf pdb outGrid

    # Load the molecule.
    set dnaMol [mol load psf $psf pdb $pdb]
    # Get the transforms for each basepair.
    set transList {}
    for {set i 0} {$i < $nBasepairs} {incr i} {
	set a [expr {$resA + $i}]
	set b [expr {$resB - $i}]
	set trans [getBasepairTransform $segA $a $segB $b $dnaMol]
	lappend transList $trans
	puts "TRANSFORM $i: $trans"
    }

    # Read the grids.
    readDx pmf $pmfGrid
    readDx local $localGrid
    set n $local(size)

    set complete 0
    for {set j 0} {$j < $n} {incr j} {
	set pos [indexToWorld local $j]
	
	# Get the value for each basepair.
	set pmfSum 0.0
	foreach trans $transList {
	    foreach {basis origin} $trans { break }

	    # Transform pos into the space of the helical DNA for this basepair.
	    set r [vecAdd [vecTransform $basis $pos] $origin]
	    
	    # Get the value of the PMF here.
	    set val [interpolatePotential pmf $r]
	    set pmfSum [expr {$pmfSum + $val}]
	}
	
	# Get the average PMF.
	lset local(data) $j [expr {$pmfSum/$nBasepairs}]

	set co [expr {int((100.0*$j)/$n)}]
	if {$co > $complete} {
	    set complete $co
	    puts "$complete percent complete"
	    writeDx local preview.dx
	}
    }

    # Write the result.
    writeDx local $outGrid
    return
}

# Get the position of an atom.
proc getPos {segName resId name {mole top}} {
    set sel [atomselect $mole "segname $segName and resid $resId and name $name"]
    set n [$sel num]
    
    if {$n < 1} {
	puts "Warning! Atom ${segName}:${resId}:${name} does not exist."
	return [list 0.0 0.0 0.0]
    } elseif {$n > 1} {
	puts "Warning! Atom ${segName}:${resId}:${name} in not unique."
    }

    set r [lindex [$sel get {x y z}] 0]
    $sel delete
    return $r
}

# Define a basis for a base.
proc getBaseNormal {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
    # Get the residue name.
    set sel [atomselect $mole $selText]
    set resName [lindex [$sel get resname] 0]
    $sel delete

    # Get the hexagon ring basis.
    if {[string equal ADE $resName] || [string equal GUA $resName]} {
	set selX0 [atomselect $mole "($selText) and name C4"]
	set selX1 [atomselect $mole "($selText) and name N1"]
	set selY0 [atomselect $mole "($selText) and name N3"]
	set selY1 [atomselect $mole "($selText) and name C5"]
    } else {
	set selX0 [atomselect $mole "($selText) and name N1"]
	set selX1 [atomselect $mole "($selText) and name C4"]
	set selY0 [atomselect $mole "($selText) and name C2"]
	set selY1 [atomselect $mole "($selText) and name C6"]
    }

    set rX0 [lindex [$selX0 get {x y z}] 0]
    set rX1 [lindex [$selX1 get {x y z}] 0]
    set rY0 [lindex [$selY0 get {x y z}] 0]
    set rY1 [lindex [$selY1 get {x y z}] 0]

    $selX0 delete
    $selX1 delete
    $selY0 delete
    $selY1 delete

    set ex [vecsub $rX1 $rX0]
    set ex [vecscale [expr 1.0/[veclength $ex]] $ex]
    set ey [vecsub $rY1 $rY0]
    set ey [vecsub $ey [vecscale [vecdot $ey $ex] $ex]]
    set ey [vecscale [expr 1.0/[veclength $ey]] $ey]
    set ez [veccross $ex $ey]

    return $ez
}

# Define a basis for a base.
proc getBasepairDirection {segA resA segB resB {mole top}} {
    set car1A [getPos $segA $resA "C1'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    return [vecUnit [vecsub $car1B $car1A]]
}

proc getBasepairBasisBest {segA resA segB resB {mole top}} {
    set zA [getBaseNormal $segA $resA $mole]
    set zB [getBaseNormal $segB $resB $mole]
    set ez [vecUnit [vecsub $zA $zB]]

    set ex [getBasepairDirection $segA $resA $segB $resB $mole]
    set ex [vecUnit [vecsub $ex [vecscale [vecdot $ex $ez] $ez]]]
    
    set ey [veccross $ez $ex]
    return [matTranspose [list $ex $ey $ez]]
}

# Get the standard position of a DNA basepair.
proc getBasepairPos {segA resA segB resB {mole top}} {
    set car1A [getPos $segA $resA "C1'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    return [vecscale 0.5 [vecadd $car1A $car1B]] 
}

extractBasepairPmf
exit
