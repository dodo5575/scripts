# This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>


set segA ADNA
set segB BDNA
# Attach the tail residues to the head residues.
set tailA 43
set tailB 20
set step 1
#Input:
set psf dna62eq.psf
set pdb dna62_splice0.pdb
#Output:
set outPdb dna62_splice1.pdb

source vector.tcl
# Standard per basepair rotation and displacement for ideal dsDNA.
set basepairRotate {{0.823022550019 -0.568008699018 5.56952340384e-07} {0.568008699021 0.823022550016 -1.59495711725e-07} {-3.6778938249e-07 4.47622281064e-07 1.0}}
set basepairDisplace {-0.471429278054 1.40150010586 3.39539747479}

if {$step < 0} {
    set tmp $tailA
    set tailA $tailB
    set tailB $tmp
    set tmp $segA
    set segA $segB
    set segB $tmp
    set basepairRotate [matTranspose $basepairRotate]
}
set headA [expr $tailA-1]
set headB [expr $tailB+1]
set attachText "(segname $segA and resid >= $tailA) or (segname $segB and resid <= $tailB)"

# Get the position of an atom.
proc getPos {segName resId name} {
    set sel [atomselect top "segname $segName and resid $resId and name $name"]
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
proc getBasepairPos {segA resA segB resB} {
    set car1A [getPos $segA $resA "C1'"]
    set car1B [getPos $segB $resB "C1'"]

    return [vecScale 0.5 [vecAdd $car1A $car1B]] 
}


# Get the standard basis of a DNA basepair.
proc getBasepairBasis {segA resA segB resB} {
    # Get the atom positions.
    set oxy5A [getPos $segA $resA "O5'"]
    set oxy3A [getPos $segA $resA "O3'"]
    set car1A [getPos $segA $resA "C1'"]
    
    set oxy5B [getPos $segB $resB "O5'"]
    set oxy3B [getPos $segB $resB "O3'"]
    set car1B [getPos $segB $resB "C1'"]

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


# Load the molecule.
mol load psf $psf pdb $pdb

# Get the positions and bases.
set headBasis [getBasepairBasis $segA $headA $segB $headB]
set headPos [getBasepairPos $segA $headA $segB $headB]
set tailBasis [getBasepairBasis $segA $tailA $segB $tailB]
set tailPos [getBasepairPos $segA $tailA $segB $tailB]

# Transform the selection to its local basis.
# That is, move tailPos=origin and tailBasis=identity.
set sel [atomselect top $attachText]
$sel moveby [vecinvert $tailPos]
$sel move [matMake4 [matTranspose $tailBasis]]

# Now attach the selection to the head.
$sel move [matMake4 [matMul $basepairRotate $headBasis]]
$sel moveby [vecAdd $headPos [vecTransform $headBasis $basepairDisplace]]
$sel delete

# Write the result.
set all [atomselect top all]
$all writepdb $outPdb
$all delete

mol delete top
exit



