# Use with: vmd -dispdev text -e generateHairpin.tcl
# Author: Jeff Comer <jcomer2@uiuc.edu>

set resA 24; # index of helix positive start
set resB 88; # index of helix negative start
set nBasepairs 7; # length of helix
# Input:
set psf dna_G-C.psf
set pdb dna_G-C.pdb
set segA ADNA
set segB BDNA
set nucleoDir dna
set nucleoSuffix "_char"
# Output:
set outName dna_G-C_fix

source vector.tcl

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

# Get the standard position of a DNA basepair.
proc getBasepairPos {segA resA segB resB {mole top}} {
    set car1A [getPos $segA $resA "C1'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    return [vecScale 0.5 [vecAdd $car1A $car1B]] 
}


# Get the standard basis of a DNA basepair.
proc getBasepairBasis {segA resA segB resB {mole top}} {
    set selText "segname $segA and resid $resA"
    set sel [atomselect $mole "segname $segA and resid $resA"]
    set resName [lindex [$sel get resname] 0]
    $sel delete

    # Get the hexagon ring basis of strand A.
    if {[string equal ADE $resName] || [string equal GUA $resName]} {
	set selX0 [atomselect $mole "segname $segA and resid $resA and name C4"]
	set selX1 [atomselect $mole "segname $segA and resid $resA and name N1"]
	set selY0 [atomselect $mole "segname $segA and resid $resA and name C4"]
	set selY1 [atomselect $mole "segname $segA and resid $resA and name C6"]
    } else {
	set selX0 [atomselect $mole "segname $segA and resid $resA and name C6"]
	set selX1 [atomselect $mole "segname $segA and resid $resA and name N3"]
	set selY0 [atomselect $mole "segname $segA and resid $resA and name C6"]
	set selY1 [atomselect $mole "segname $segA and resid $resA and name C4"]
    }
    if {[$selX0 num] < 1 || [$selX1 num] < 1 || [$selY0 num] < 1 || [$selY1 num] < 1} {
	puts stderr "A necessary atom for $segA:$resA does not exist."
	return [matIdentity]
    }
    set dxA [vecsub [lindex [$selX1 get {x y z}] 0] [lindex [$selX0 get {x y z}] 0]]
    set dyA [vecsub [lindex [$selY1 get {x y z}] 0] [lindex [$selY0 get {x y z}] 0]]
    $selX0 delete
    $selX1 delete
    $selY0 delete
    $selY1 delete

    # Get the hexagon ring basis of strand B.
    # dxB ~ dxA and dyB ~ dyA for paired bases.
    if {[string equal ADE $resName] || [string equal GUA $resName]} {
	set selX0 [atomselect $mole "segname $segB and resid $resB and name N3"]
	set selX1 [atomselect $mole "segname $segB and resid $resB and name C6"]
	set selY0 [atomselect $mole "segname $segB and resid $resB and name N3"]
	set selY1 [atomselect $mole "segname $segB and resid $resB and name C5"]
    } else {
	set selX0 [atomselect $mole "segname $segB and resid $resB and name N1"]
	set selX1 [atomselect $mole "segname $segB and resid $resB and name C4"]
	set selY0 [atomselect $mole "segname $segB and resid $resB and name N1"]
	set selY1 [atomselect $mole "segname $segB and resid $resB and name N3"]
    }
    if {[$selX0 num] < 1 || [$selX1 num] < 1 || [$selY0 num] < 1 || [$selY1 num] < 1} {
	puts stderr "A necessary atom for $segB:$resB does not exist."
	return [matIdentity]
    }
    set dxB [vecsub [lindex [$selX1 get {x y z}] 0] [lindex [$selX0 get {x y z}] 0]]
    set dyB [vecsub [lindex [$selY1 get {x y z}] 0] [lindex [$selY0 get {x y z}] 0]]
    $selX0 delete
    $selX1 delete
    $selY0 delete
    $selY1 delete

    # Use the mean of the base bases.
    set ex [vecscale 0.5 [vecadd $dxA $dxB]]
    #set ex $dxB
    set ex [vecscale [expr 1.0/[veclength $ex]] $ex]
    set ey [vecscale 0.5 [vecadd $dyA $dyB]]
    #set ey $dyB
    set ey [vecsub $ey [vecscale [vecdot $ey $ex] $ex]]
    set ey [vecscale [expr 1.0/[veclength $ey]] $ey]
    set ez [veccross $ex $ey]

    return [matTranspose [list $ex $ey $ez]]
}

proc getBasepairBasis1 {segA resA segB resB {mole top}} {
    # Get the atom positions.
    set oxy5A [getPos $segA $resA "O5'" $mole]
    set oxy3A [getPos $segA $resA "O3'" $mole]
    set car1A [getPos $segA $resA "C1'" $mole]
    
    set oxy5B [getPos $segB $resB "O5'" $mole]
    set oxy3B [getPos $segB $resB "O3'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    # The first basis vector is along O5' -> O3'.
    set oxy5 [vecScale 0.5 [vecAdd $oxy5A $oxy5B]] 
    set oxy3 [vecScale 0.5 [vecAdd $oxy3A $oxy3B]]
    set a [vecSub $oxy3 $oxy5]
    set a [vecScale [expr 1.0/[vecLength $a]] $a]
    
    # The second basis vector tends along A:C1' -> B:C1'.
    # Remove the component along $a.
    set b [vecSub $car1B  $car1A]
    set b [vecSub $b [vecScale [vecDot $a $b] $a]]
    set b [vecScale [expr 1.0/[vecLength $b]] $b]
    
    # The third basis vector is $a cross $b.
    set c [vecCross $a $b]

    return [matTranspose [list $a $b $c]]
}

proc fitBasepairs {segA resA segB resB atMol gcMol} {
    set any 1;

    set sel [atomselect top "segname $segA and resid $resA"]
    set resName [lindex [$sel get resname] 0]

    if {[string equal ADE $resName]} {
	set pairMol $atMol
    } elseif {[string equal THY $resName]} {
	set pairMol $atMol
	set tmp $segA
	set segA $segB
	set segB $tmp
	set tmp $resA
	set resA $resB
	set resB $tmp
    } elseif {[string equal GUA $resName]} {
	set pairMol $gcMol
    } else {
	set pairMol $gcMol
	set tmp $segA
	set segA $segB
	set segB $tmp
	set tmp $resA
	set resA $resB
	set resB $tmp
    }

    set basis [getBasepairBasis $segA $resA $segB $resB]
    set pos [getBasepairPos $segA $resA $segB $resB]
    set basis0 [getBasepairBasis ADNA $any BDNA $any $pairMol]
    set pos0 [getBasepairPos ADNA $any BDNA $any $pairMol]

    # Transform strand A.
    set selA0 [atomselect $pairMol "segname ADNA and resid $any"]
    foreach name [lsort -unique [$selA0 get name]] {
	set s0 [atomselect $pairMol "segname ADNA and resid $any and name $name"]
	set s [atomselect top "segname $segA and resid $resA and name $name"]
	if {[$s num] < 1} {
	    puts stderr "Warning! $segA:$resA:$name does not exist."
	    continue
	}

	set r [lindex [$s0 get {x y z}] 0]
	set r [vecTransform [matTranspose $basis0] [vecsub $r $pos0]]
	set r [vecadd [vecTransform $basis $r] $pos]
	$s set {x y z} [list $r]
	$s0 delete
	$s delete
    }
    
    # Transform strand B.
    set selB0 [atomselect $pairMol "segname BDNA and resid $any"]
    foreach name [lsort -unique [$selB0 get name]] {
	set s0 [atomselect $pairMol "segname BDNA and resid $any and name $name"]
	set s  [atomselect top "segname $segB and resid $resB and name $name"]
	if {[$s num] < 1} {
	    puts stderr "Warning! $segB:$resB:$name does not exist."
	    continue
	}

	set r [lindex [$s0 get {x y z}] 0]
	set r [vecTransform [matTranspose $basis0] [vecsub $r $pos0]]
	set r [vecadd [vecTransform $basis $r] $pos]
	$s set {x y z} [list $r]
	$s0 delete
	$s delete
    }
    
    return
}


# Load the basepair templates.
set atMol [mol load pdb $nucleoDir/at${nucleoSuffix}.pdb]
set gcMol [mol load pdb $nucleoDir/gc${nucleoSuffix}.pdb]

# Load the molecule.
mol load psf $psf pdb $pdb

# Make basepairs match up.
for {set i 0} {$i < $nBasepairs} {incr i} {
    set a [expr $resA + $i]
    set b [expr $resB - $i]

    puts "$a $b"
    fitBasepairs $segA $a $segB $b $atMol $gcMol
}

set all [atomselect top all]
$all writepsf $outName.psf
$all writepdb $outName.pdb

mol delete $gcMol
mol delete $atMol
mol delete top
exit
