# Use with: vmd -dispdev text -e templateDna.tcl
# Author: Jeff Comer <jcomer2@uiuc.edu>

# Input:
set templatePsf coated_trap_charmm.psf
set templateCoord coated_trap_charmm.pdb
set mapPsf coated_trap_charmm.psf
set mapCoord  trap_coated_4V0.restart.coor
set selText "segname ADNA"
set nucleoDir "dna"
set nucleoSuffix "_char"
# Output:
set pdb dna_fix.pdb

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
proc mapCoordinates {templateMol mapMol seg res nameList} {
    # Handle nameList == all.
    if {[string equal $nameList all]} {
	set s [atomselect $mapMol "segname $seg and resid $res"]
	set nameList [$s get name]
	$s delete
    }
    set missList {}

    foreach nam $nameList {
	set sel0 [atomselect $templateMol "segname $seg and resid $res and name $nam"]
	set sel1 [atomselect $mapMol "segname $seg and resid $res and name $nam"]

	if {[$sel0 num] < 1} {
	    puts "WARNING! Cannot find $seg:$res:$nam in molecule $templateMol"
	    lappend missList [list $seg $res $nam]
	    continue
	} elseif {[$sel0 num] > 1} {
	    puts "WARNING! $seg:$res:$nam in molecule $templateMol is not unique"
	    lappend missList [list $seg $res $nam]
	    continue
	}

	if {[$sel1 num] < 1} {
	    puts "WARNING! Cannot find $seg:$res:$nam in molecule $mapMol"
	    lappend missList [list $seg $res $nam]
	    continue
	} elseif {[$sel1 num] > 1} {
	    puts "WARNING! $seg:$res:$nam in molecule $mapMol is not unique"
	    lappend missList [list $seg $res $nam]
	    continue
	}

	set pos [$sel0 get {x y z}]
	$sel1 set {x y z} $pos
    }

    return $missList
}

# Map coordinates from one molecule to another.
proc mapCoordinatesTemp {templateMol mapMol seg res nameList} {
    # Handle nameList == all.
    if {[string equal $nameList all]} {
	set s [atomselect $mapMol "segname $seg and resid $res"]
	set nameList [$s get name]
	$s delete
    }
    set missList {}

    foreach nam $nameList {
	set sel0 [atomselect $templateMol "name $nam"]
	set sel1 [atomselect $mapMol "segname $seg and resid $res and name $nam"]

	if {[$sel0 num] < 1} {
	    puts "WARNING! Cannot find $nam in molecule $templateMol"
	    lappend missList [list $seg $res $nam]
	    continue
	} elseif {[$sel0 num] > 1} {
	    puts "WARNING! $nam in molecule $templateMol is not unique"
	    lappend missList [list $seg $res $nam]
	    continue
	}

	if {[$sel1 num] < 1} {
	    puts "WARNING! Cannot find $seg:$res:$nam in molecule $mapMol"
	    lappend missList [list $seg $res $nam]
	    continue
	} elseif {[$sel1 num] > 1} {
	    puts "WARNING! $seg:$res:$nam in molecule $mapMol is not unique"
	    lappend missList [list $seg $res $nam]
	    continue
	}

	set pos [$sel0 get {x y z}]
	$sel1 set {x y z} $pos
    }

    return $missList
}

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
set templateMol [mol load psf $templatePsf]
mol addfile $templateCoord
puts "Loaded the template molecule."

set mapMol [mol load psf $mapPsf]
mol addfile $mapCoord
set mapSel [atomselect $mapMol $selText]
foreach zero {0} {set resList [lsort -unique [$mapSel get {segname resid resname}]]}
puts "Loaded the molecule to be mapped."

puts "Mapping the nucleotides..."
set missList {}
foreach resSet $resList {
    foreach {segName resId resName} $resSet {break}

    # Check the residue name.
    set good 0
    foreach res $resNameList {
	if {[string equal $res $resName]} {
	    set good 1
	    break
	}
    }
    if {!$good} {
	puts "Warning! Unrecognized residue: $resName."
	continue
    }
    puts "Mapping $segName:$resName:$resId"
    
    # Get the template basis and position.
    set templateBasis [getBaseBasis $segName $resId $templateMol]
    set templatePos [getBasePos $segName $resId $templateMol]

    # Map the coordinates from the ideal nucleotide.
    set missList [concat $missList [mapCoordinatesTemp $loadMol($resName) $mapMol $segName $resId all]]

    # Shift them into the template basis.
    set s [atomselect $mapMol "segname $segName and resid $resId"]
    $s moveby [vecInvert $loadPos($resName)]
    $s move [matMake4 [matTranspose $loadBasis($resName)]]
    $s move [matMake4 $templateBasis]
    $s moveby $templatePos
}

# Map the missed atoms to the template.
puts "\nFixing any missed atoms if they exist in the template..."
foreach atom $missList {
    foreach {segName resId name} $atom {break}
    
    puts "Mapping atom $segName:$resId:$name"
    mapCoordinates $templateMol $mapMol $segName $resId [list $name]
}

set all [atomselect $mapMol all]
$all writepdb $pdb
$all delete

$mapSel delete
foreach res $resNameList {
    mol delete $loadMol($res)
}
mol delete $mapMol
mol delete $templateMol

exit
