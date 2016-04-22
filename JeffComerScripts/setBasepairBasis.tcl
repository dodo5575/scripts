# Author: jcomer2@illinois.edu

set segA ADNA
set resA 1
set segB BDNA
set resB 60
set nBasepairs 60
set selText nucleic
# Input:
set psf trap_periodic.psf
set pdb period_trap_0Va6.pdb
# Output:
set outName trap_basis

source $env(HOME)/scripts/vector.tcl

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
    set sel [atomselect top $selText]
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

mol load psf $psf pdb $pdb

# Make a list of the residues.
set resAList {}
set segAList {}
set resBList {}
set segBList {}
for {set i 0} {$i < $nBasepairs} {incr i} {
    set a [expr {$resA + $i}]
    set b [expr {$resB - $i}]
    lappend segAList $segA
    lappend resAList $a
    lappend segBList $segB
    lappend resBList $b

    set selA [atomselect top "segname $segA and resid $a"]
    set selB [atomselect top "segname $segB and resid $b"]
    if {[$selA num] < 3} {
	puts "ERROR! Residue $segA:$a does not exist."
	exit
    }
    if {[$selB num] < 3} {
	puts "ERROR! Residue $segB:$b does not exist."
	exit
    }

    puts "BASEPAIR $segA:$a $segB:$b"
    $selA delete
    $selB delete
}

foreach segA $segAList resA $resAList segB $segBList resB $resBList {
    set sel [atomselect top "(segname $segA and resid $resA) or (segname $segB and resid $resB)"]
    set pos [getBasepairPos $segA $resA $segB $resB]
    set rot [getBasepairBasisBest $segA $resA $segB $resB]
    
    $sel moveby [vecinvert $pos]
    $sel move [matMake4 [matTranspose $rot]]
    #$sel moveby $pos
    $sel delete
}

set all [atomselect top all]
$all writepsf $outName.psf
$all writepdb $outName.pdb

exit
