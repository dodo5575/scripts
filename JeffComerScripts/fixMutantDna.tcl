# Use with: vmd -dispdev text -e templateDna.tcl
# Author: Jeff Comer <jcomer2@uiuc.edu>
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

# Get the standard position of DNA nucleotide.
proc getNucleotidePos {seg res {mole top}} {
    set rC5P [getPos $seg $res "C5'" $mole]
    set rC3P [getPos $seg $res "C3'" $mole]

    return [vecScale 0.5 [vecAdd $rC3P $rC5P]] 
}

# Define a basis for the nucleotide.
proc getNucleotideBasis {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
    set selO5P [atomselect $mole "($selText) and name O5'"]
    set selO3P [atomselect $mole "($selText) and name O3'"]
    set selC1P [atomselect $mole "($selText) and name C1'"]

    set resName [lindex [$selO5P get resname] 0]
    if {[string equal ADE $resName] || [string equal GUA $resName]} {
	set selN [atomselect $mole "($selText) and name N9"]
    } else {
	set selN [atomselect $mole "($selText) and name N1"]
    }

    set rO5P [lindex [$selO5P get {x y z}] 0]
    set rO3P [lindex [$selO3P get {x y z}] 0]
    set rC1P [lindex [$selC1P get {x y z}] 0]
    set rN [lindex [$selN get {x y z}] 0]
    
    $selO5P delete
    $selO3P delete
    $selC1P delete
    $selN delete

    set ex [vecsub $rO3P $rO5P]
    set ex [vecscale [expr 1.0/[veclength $ex]] $ex]
    set ey [vecsub $rC1P $rO5P]
    set ey [vecsub $ey [vecscale [vecdot $ey $ex] $ex]]
    set ey [vecscale [expr 1.0/[veclength $ey]] $ey]
    set ez [veccross $ex $ey]

    return [matTranspose [list $ex $ey $ez]]
}
# Get the standard position of DNA base.
proc getBasePos {seg res {mole top}} {
    return [getPos $seg $res "C1'" $mole]
}

# Define a basis for the base.
proc getBaseBasis {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
    set sel [atomselect $mole $selText]
    set resName [lindex [$sel get resname] 0]
    $sel delete

    # Get the hexagon ring basis.
    if {[string equal ADE $resName] || [string equal GUA $resName]} {
	set selX0 [atomselect $mole "($selText) and name C4"]
	set selX1 [atomselect $mole "($selText) and name N1"]
	set selY0 [atomselect $mole "($selText) and name C4"]
	set selY1 [atomselect $mole "($selText) and name C6"]
    } else {
	set selX0 [atomselect $mole "($selText) and name C6"]
	set selX1 [atomselect $mole "($selText) and name N3"]
	set selY0 [atomselect $mole "($selText) and name C6"]
	set selY1 [atomselect $mole "($selText) and name C4"]
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

# Map coordinates from one molecule to another.
proc mapCoordinates {tempMol tempSeg tempRes targMol targSeg targRes nameList} {
    # Handle nameList == all.
    if {[string equal $nameList all]} {
	set s [atomselect $targMol "segname $targSeg and resid $targRes"]
	set nameList [$s get name]
	$s delete
    }
    set missList {}

    foreach nam $nameList {
	set sel0 [atomselect $tempMol "segname $tempSeg and resid $tempRes and name $nam"]
	set sel1 [atomselect $targMol "segname $targSeg and resid $targRes and name $nam"]

	if {[$sel0 num] < 1} {
	    #puts "WARNING! Cannot find $tempSeg:$tempRes:$nam in molecule $tempMol"
	    lappend missList [list $tempSeg $tempRes $nam]
	    continue
	} elseif {[$sel0 num] > 1} {
	    puts "WARNING! $tempSeg:$tempRes:$nam in molecule $tempMol is not unique"
	    lappend missList [list $tempSeg $tempRes $nam]
	    continue
	}

	if {[$sel1 num] < 1} {
	    puts "WARNING! Cannot find $targSeg:$targRes:$nam in molecule $targMol"
	    lappend missList [list $targSeg $targRes $nam]
	    continue
	} elseif {[$sel1 num] > 1} {
	    puts "WARNING! $targSeg:$targRes:$nam in molecule $targMol is not unique"
	     lappend missList [list $targSeg $targRes $nam]
	    continue
	}

	set pos [$sel0 get {x y z}]
	$sel1 set {x y z} $pos
    }

    return $missList
}

# Map coordinates from one molecule to another.
proc mapCoordinatesTemp {tempMol targMol seg res nameList} {
    # Handle nameList == all.
    if {[string equal $nameList all]} {
	set s [atomselect $targMol "segname $seg and resid $res"]
	set nameList [$s get name]
	$s delete
    }
    set missList {}

    foreach nam $nameList {
	set sel0 [atomselect $tempMol "name $nam"]
	set sel1 [atomselect $targMol "segname $seg and resid $res and name $nam"]

	if {[$sel0 num] < 1} {
	    puts "WARNING! Cannot find $nam in molecule $tempMol"
	    lappend missList [list $seg $res $nam]
	    continue
	} elseif {[$sel0 num] > 1} {
	    puts "WARNING! $nam in molecule $tempMol is not unique"
	    lappend missList [list $seg $res $nam]
	    continue
	}

	if {[$sel1 num] < 1} {
	    puts "WARNING! Cannot find $seg:$res:$nam in molecule $targMol"
	    lappend missList [list $seg $res $nam]
	    continue
	} elseif {[$sel1 num] > 1} {
	    puts "WARNING! $seg:$res:$nam in molecule $targMol is not unique"
	    lappend missList [list $seg $res $nam]
	    continue
	}

	set pos [$sel0 get {x y z}]
	$sel1 set {x y z} $pos
    }

    return $missList
}

proc fixMutantDna {tempPdb tempText targPdb targText outPdb} {
    set nucleoDir "dna"
    set nucleoSuffix "_char"
    set backboneAtoms "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P H3T H5T"

    # Load the ideal nucleotide templates.
    set resNameList {ADE CYT GUA THY}
    foreach res $resNameList {
	set res1 [string tolower $res]

	set loadMol($res) [mol load pdb $nucleoDir/${res1}${nucleoSuffix}.pdb]
	set sel [atomselect top all]
	set segName [lindex [$sel get segname] 0]
	set resId [lindex [$sel get resid] 0]
	set loadBasis($res) [getBaseBasis $segName $resId]
	set loadPos($res) [getBasePos $segName $resId]
	set loadResName($res) [lindex [$sel get resname] 0]
	set loadName($res) [$sel get name]
	$sel delete
    }
    puts "Loaded ideal nucleotide templates."


    # Load the molecule.
    set tempMol [mol load pdb $tempPdb]
    set tempSel [atomselect $tempMol $tempText]
    set tempResList [lsort -unique [$tempSel get {segname resid resname}]]
    set tempResList [lsort -index 1 -integer $tempResList]
    puts "Loaded the template molecule $tempPdb."
    puts "Selected [llength $tempResList] residues."
    puts "First residue: [lindex $tempResList 0]"
    
    set targMol [mol load pdb $targPdb]
    set targSel [atomselect $targMol $targText]
    set targResList [lsort -unique [$targSel get {segname resid resname}]]
    set targResList [lsort -index 1 -integer $targResList]
    puts "Loaded the target molecule $targPdb."
    puts "Selected [llength $targResList] residues."
    puts "First residue: [lindex $targResList 0]"
    set nRes [llength $targResList]
    
    if {$nRes > [llength $tempResList]} {
	puts "Error! There are more residues in the target than in the template!"
	return -1
    }
    
    puts "\nMapping the bases and backbones..."
    for {set i 0} {$i < $nRes} {incr i} {
	foreach {tempSeg tempRes tempResName} [lindex $tempResList $i] { break }
	foreach {targSeg targRes targResName} [lindex $targResList $i] { break }

	# Check the residue name.
	set good 0
	foreach res $resNameList {
	    if {[string equal $res $targResName]} {
		set good 1
		break
	    }
	}
	if {!$good} {
	    puts "Warning! Unrecognized residue: $targResName."
	    continue
	}
	
	if {[string equal $targResName $tempResName]} {
	    # Map the atoms exactly if the residue names are the same.
	    puts "exact fit: target $targSeg:$targRes:$targResName -> template $tempSeg:$tempRes:$tempResName"
	    # Map all the atoms possible.
	    mapCoordinates $tempMol $tempSeg $tempRes $targMol $targSeg $targRes all

	} else {
	    # Map based on a fit otherwise.
	    puts "base fit: target $targSeg:$targRes:$targResName -> template $tempSeg:$tempRes:$tempResName"
	    # Map the coordinates from the ideal nucleotide.
	    mapCoordinatesTemp $loadMol($targResName) $targMol $targSeg $targRes all

	    # Get the template basis and position.
	    set tempBasis [getBaseBasis $tempSeg $tempRes $tempMol]
	    set tempPos [getBasePos $tempSeg $tempRes $tempMol]
	    
	    # Shift into the template basis.
	    set s [atomselect $targMol "segname $targSeg and resid $targRes"]
	    $s moveby [vecInvert $loadPos($targResName)]
	    $s move [matMake4 [matTranspose $loadBasis($targResName)]]
	    $s move [matMake4 $tempBasis]
	    $s moveby $tempPos
	    
	    # Map the backbone.
	    mapCoordinates $tempMol $tempSeg $tempRes $targMol $targSeg $targRes $backboneAtoms
	}
    }

    set all [atomselect $targMol all]
    $all writepdb $outPdb
    $all delete
    puts "\nWrote $outPdb."

    $targSel delete
    foreach res $resNameList {
	mol delete $loadMol($res)
    }
    mol delete $targMol
    mol delete $tempMol

    return 0
}
