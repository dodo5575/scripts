# Use with: vmd -dispdev text -e generateHairpin.tcl
# Author: Jeff Comer <jcomer2@uiuc.edu>

set resA 24; # index of helix positive start
set resB 89; # index of helix negative start
#set resA 27; # index of helix positive start
#set resB 84; # index of helix negative start
set nBasepairs 9; # length of helix
set nTries 5000
# Input:
set templatePsf cut_dna_patched0.psf
set templatePdb cut_dna_patched0.pdb
set psf dna_A-T_fix.psf
set pdb dna_A-T_map.pdb
set segA ADNA
set segB BDNA
# Output:
set outName dna_A-T_better

set backboneAtoms "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"
source vector.tcl

# Define the hydrogen-bond-forming atoms.
set hydrogenBondLight {{ADE N1 THY H3 1.81} {ADE H61 THY O4 1.87} {GUA H21 CYT O2 1.80} {GUA H1 CYT N3 1.84} {GUA O6 CYT H41 1.80}}
set hydrogenBond {{ADE N1 THY N3 2.90} {ADE N6 THY O4 2.88} {GUA N2 CYT O2 2.81} {GUA N1 CYT N3 2.87} {GUA O6 CYT N4 2.81}}
# Make the hydrogen bond list symmetric by adding reverses.
foreach b $hydrogenBond {
    lappend hydrogenBond [list [lindex $b 2] [lindex $b 3] [lindex $b 0] [lindex $b 1] [lindex $b 4]]}
# Make list of the first resname.
set hBondResName {}
foreach b $hydrogenBond {
    lappend hBondResName [lindex $b 0]
}

proc mapCoordinates {templateMol mapMol seg res nameList} {
    foreach nam $nameList {
	set sel0 [atomselect $templateMol "segname $seg and resid $res and name $nam"]
	set sel1 [atomselect $mapMol "segname $seg and resid $res and name $nam"]
	
	if {[$sel0 num] < 1} {
	    puts "WARNING! Cannot find $seg:$res:$nam in molecule $templateMol"
	    return
	} elseif {[$sel0 num] > 1} {
	     puts "WARNING! $seg:$res:$nam in molecule $templateMol is not unique"

	}

	if {[$sel1 num] < 1} {
	    puts "WARNING! Cannot find $seg:$res:$nam in molecule $mapMol"
	    return -1
	} elseif {[$sel1 num] > 1} {
	     puts "WARNING! $seg:$res:$nam in molecule $mapMol is not unique"
	}

	set pos [$sel0 get {x y z}]
	$sel1 set {x y z} $pos
    }

    return 0
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

# Get the standard position of DNA base.
proc getBasePos {seg res {mole top}} {
    return [getPos $seg $res "C1'" $mole]
}

# Define a basis for a base.
proc getBaseBasis {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
 
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

    return [matTranspose [list $ex $ey $ez]]
}

# Get the standard position of a DNA basepair.
proc getBasepairPos {segA resA segB resB {mole top}} {
    set car1A [getPos $segA $resA "C1'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    return [vecScale 0.5 [vecAdd $car1A $car1B]] 
}

# Get the standard basis of a DNA basepair.
proc getBasepairBasis {segA resA segB resB {mole top}} {
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

proc getHydrogenBonds {segA resA segB resB {mole top}} {
    global hydrogenBond hBondResName
    # Get the residue name.
    set sel [atomselect $mole "segname $segA and resid $resA"]
    set resName [lindex [$sel get resname] 0]
    $sel delete

    set hbA {}
    set hbB {}
    set hbDist {}
    set indList [lsearch -all $hBondResName $resName]
    foreach ind $indList {
	foreach {baseA nameA baseB nameB dist} [lindex $hydrogenBond $ind] {break}

	set selA [atomselect $mole "segname $segA and resid $resA and name $nameA"]
	set selB [atomselect $mole "segname $segB and resid $resB and name $nameB"]
	if {[$selA num] != 1 && [$selB num] != 1}  {
	    puts stderr "Warning! Cannot bond $baseA to $baseB."
	    puts stderr "Atoms $resA:$nameA and $resB:$nameB were not found."
	} else {
	    set iA [lindex [$selA get index] 0]
	    set iB [lindex [$selB get index] 0]
	    lappend hbA $iA
	    lappend hbB $iB
	    lappend hbDist $dist
	    incr nBonds
	}
	$selA delete
	$selB delete
    }

    return [list $hbA $hbB $hbDist]
}

proc hydrogenBondEnergy {hb} {
    foreach {hbA hbB hbDist} $hb {break}

    set e 0.0
    foreach a $hbA b $hbB d0 $hbDist {
	set d [measure bond [list $a $b]]
	#set e [expr $e + ($d-$d0)*($d-$d0)]
	set e [expr $e + $d*$d]
    }
    return $e
}

# Load the molecule.
set templateMol [mol load psf $templatePsf pdb $templatePdb]
set mapMol [mol load psf $psf pdb $pdb]

for {set i 0} {$i < $nBasepairs} {incr i} {
    set a [expr $resA + $i]
    set b [expr $resB - $i]

    puts "$a $b"

    # Get the B strand bases.
    set mapPos [getBasePos $segB $b $mapMol]
    set mapBasis [getBaseBasis $segB $b $mapMol]
    set tempPos [getBasePos $segB $b $templateMol]
    set tempBasis [getBaseBasis $segB $b $templateMol]

    # Map the B strand bases.
    set selB [atomselect $mapMol "segname $segB and resid $b"]
    $selB moveby [vecInvert $mapPos]
    $selB move [matMake4 [matTranspose $mapBasis]]
    $selB move [matMake4 $tempBasis]
    $selB moveby $tempPos

    # Get the origin.
    set originSel [atomselect $mapMol "segname $segA and resid $a and name C1'"]
    set origin [lindex [$originSel get {x y z}] 0]
    $originSel delete
    
    # Get the hydrogen bonds.
    set hb [getHydrogenBonds $segA $a $segB $b $mapMol]

    # Try random rotations of the A strand to get the best one.
    set sel [atomselect $mapMol "segname $segA and resid $a"]
    set coord [$sel get {x y z}]
    set bestRotate [matIdentity]
    set bestEnergy [hydrogenBondEnergy $hb]
    for {set j 0} {$j < $nTries} {incr j} {
	$sel set {x y z} $coord

	# Make the rotation.
	set rotate [matRandomRot]
	$sel moveby [vecInvert $origin]
	$sel move [matMake4 $rotate]
	$sel moveby $origin

	# Compare the energies.
	set energy [hydrogenBondEnergy $hb]
	if {$energy < $bestEnergy} {
	    set bestEnergy $energy
	    set bestRotate $rotate
	    puts $energy
	}
    }
    $sel set {x y z} $coord
    $sel moveby [vecInvert $origin]
    $sel move [matMake4 $bestRotate]
    $sel moveby $origin
    $sel delete

    # Map the backbone exactly.
    mapCoordinates $templateMol $mapMol $segA $a $backboneAtoms
    mapCoordinates $templateMol $mapMol $segB $b $backboneAtoms
}


set all [atomselect $mapMol all]
$all writepsf $outName.psf
$all writepdb $outName.pdb

mol delete $templateMol
mol delete $mapMol
exit
