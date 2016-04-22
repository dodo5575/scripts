# Match one selection with another.
# Use with: vmd -dispdev text -e fit.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set tempSegA ADNA
set tempResA 1
set tempSegB BDNA
set tempResB 111
set fitSegA ADN0
set fitResA 111
set fitSegB BDN0
set fitResB 1
set transText "all"
#Input:
set tempPsf dna.psf
set tempPdb dna_along_x.pdb
set fitPsf longDNA_charmm0.psf
set fitPdb longDNA_charmm0.pdb
#Output:
set finalPdb longDNA_charmm0_fit.pdb

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

set tempMol [mol load psf $tempPsf pdb $tempPdb]
set fitMol [mol load psf $fitPsf pdb $fitPdb]

# Get the positions and bases
set tempBasis [getBasepairBasis $tempSegA $tempResA $tempSegB $tempResB $tempMol]
set tempPos [getBasepairPos $tempSegA $tempResA $tempSegB $tempResB $tempMol]
set fitBasis [getBasepairBasis $fitSegA $fitResA $fitSegB $fitResB $fitMol]
set fitPos [getBasepairPos $fitSegA $fitResA $fitSegB $fitResB $fitMol]

# Transform the selection to its local basis.
# That is, move tailPos=origin and tailBasis=identity.
set transSel [atomselect $fitMol $transText]
$transSel moveby [vecInvert $fitPos]
$transSel move [matMake4 [matTranspose $fitBasis]]

# Now attach the selection to the head.
$transSel move [matMake4 $tempBasis]
$transSel moveby $tempPos
$transSel delete

# Write the result.
set all [atomselect $fitMol all]
$all writepdb $finalPdb
$all delete

mol delete $tempMol
mol delete $fitMol
exit
