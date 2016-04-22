# Use with: vmd -dispdev text -e generateHairpin.tcl
# Author: Jeff Comer <jcomer2@uiuc.edu>

# Define the sequence.
set seq G21

# Input:
set templatePsf cut_dna.psf
set templatePdb cut_dna.pdb
set templateSeg ADNA
set topFile top_all27_prot_na.inp
set nucleoDir dna
set nucleoSuffix "_char"
# Output:
set psf dna_${seq}.psf
set pdb dna_${seq}.pdb
set outSeg $templateSeg

source vector.tcl

# Expand a DNA sequence given by nucleotide letters and repeat counts
# to just nucleotide letters.
# Example: GA4TC2 -> GAAAATCC
proc expandSequence {seq} {
    set s ""
    set n [string length $seq]

    for {set i 0} {$i < $n} {incr i} {
	set c [string index $seq $i]
	set c1 [string index $seq [expr $i+1]]
	if {[string equal $c " "]} {continue}

	# Is it the symbol nucleotide?
	if {[regexp "\[ACGTUN\]" $c]} {
	    if {$i == $n-1 || ![string is integer $c1]} {
		set count 1
	    } else {
		# Collect the count.
		set num ""
		for {set j [expr $i+1]} {$j < $n} {incr j} {
		    set numeral [string index $seq $j] 
		    if {[string is integer $numeral]} {
			set num ${num}${numeral}
		    } else {
			break
		    }
		}
		set count $num
	    }

	    # Add the 
	    for {set j 0} {$j < $count} {incr j} {
		set s ${s}${c}
	    }
	}
    }
    return $s
}

# Map coordinates from one molecule to another.
proc mapCoordinates {templateMol mapMol seg res nameList} {
    # Handle nameList == all.
    if {[string equal $nameList all]} {
	set s [atomselect $mapMol "segname $seg and resid $res"]
	set nameList [$s get name]
	$s delete
    }

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

# Return a list with atom positions.
proc extractPdbCoords {pdbFile} {
    set r {}
    
    # Get the coordinates from the pdb file.
    set in [open $pdbFile r]
    foreach line [split [read $in] \n] {
	if {[string equal [string range $line 0 3] "ATOM"]} {
	    set x [string trim [string range $line 30 37]]
	    set y [string trim [string range $line 38 45]]
	    set z [string trim [string range $line 46 53]]
	    
	    lappend r [list $x $y $z]
	}
    }
    close $in
    return $r
}

# Extract all atom records from a pdb file.
proc extractPdbRecords {pdbFile} {
    set in [open $pdbFile r]
    
    set pdbLine {}
    foreach line [split [read $in] \n] {
	if {[string equal [string range $line 0 3] "ATOM"]} {
	    lappend pdbLine $line
	}
    }
    close $in	
    
    return $pdbLine
}

# Construct a pdb line from a template line, index, resId, and coordinates.
proc makePdbLine {template index segName resId r} {
    foreach {x y z} $r {break}
    set record "ATOM  "
    set si [string range [format "     %5i " $index] end-5 end]
    set temp0 [string range $template 12 21]
    set resId [string range "    $resId"  end-3 end]
    set temp1 [string range $template  26 29]
    set sx [string range [format "       %8.3f" $x] end-7 end]
    set sy [string range [format "       %8.3f" $y] end-7 end]
    set sz [string range [format "       %8.3f" $z] end-7 end]
    set temp2 [string range $template 54 71]
    set segName [string range "$segName    "  0 3]
    set tempEnd [string range $template 76 end]

    # Construct the pdb line.
    return "${record}${si}${temp0}${resId}${temp1}${sx}${sy}${sz}${temp2}${segName}${tempEnd}"
}

# Construct a pdb line from everything.
proc makePdbLineFull {index segName resId name resName r} {
    set template "ATOM     42  HN  ASN X   4     -41.083  17.391  50.684  0.00  0.00      P1    "

    foreach {x y z} $r {break}
    set record "ATOM  "
    set si [string range [format "     %5i " $index] end-5 end]
    if {[string length $name] < 4} {
	set name [string range " $name    " 0 3]
    } else {
	set name [string range $name 0 3]
    }
    set resName [string range " $resName    " 0 3]
    set temp0 " [string index $segName 0]"
    set resId [string range "    $resId"  end-3 end]
    set temp1 [string range $template  26 29]
    set sx [string range [format "       %8.3f" $x] end-7 end]
    set sy [string range [format "       %8.3f" $y] end-7 end]
    set sz [string range [format "       %8.3f" $z] end-7 end]
    set temp2 [string range $template 54 71]
    set segName [string range "$segName    "  0 3]
    set tempEnd [string range $template 76 end]

    # Construct the pdb line.
    return "${record}${si}${name}${resName}${temp0}${resId}${temp1}${sx}${sy}${sz}${temp2}${segName}${tempEnd}"
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

# Get the standard position of a DNA basepair.
proc getBasepairPos {segA resA segB resB {mole top}} {
    set car1A [getPos $segA $resA "C1'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    return [vecScale 0.5 [vecAdd $car1A $car1B]] 
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

# Define a basis for a basepair.
proc getBasepairBasisBest {segA resA segB resB {mole top}} {
    set zA [getBaseNormal $segA $resA $mole]
    set zB [getBaseNormal $segB $resB $mole]
    set ez [vecUnit [vecsub $zA $zB]]

    set ex [getBasepairDirection $segA $resA $segB $resB $mole]
    set ex [vecUnit [vecsub $ex [vecscale [vecdot $ex $ez] $ez]]]
    
    set ey [veccross $ez $ex]
    return [matTranspose [list $ex $ey $ez]]
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

# Get the sequence.
set sequence [expandSequence $seq]
set nSequence [string length $sequence]

# Load the information for the nucleotides.
foreach res {ade cyt gua thy} {
    set loadPos($res) [extractPdbCoords $nucleoDir/${res}${nucleoSuffix}.pdb]
    set loadRec($res) [extractPdbRecords $nucleoDir/${res}${nucleoSuffix}.pdb]

    mol load pdb $nucleoDir/${res}${nucleoSuffix}.pdb
    set sel [atomselect top all]
    set segName [lindex [$sel get segname] 0]
    set resId [lindex [$sel get resid] 0]
    set loadRot($res) [getBaseBasis $segName $resId]
    set loadCen($res) [getBasePos $segName $resId]
    set loadResName($res) [lindex [$sel get resname] 0]
    $sel delete
    mol delete top
}

# Load the template molecule.
mol load psf $templatePsf pdb $templatePdb
set out [open $outSeg.pdb w]
set n 1
set nSeq 0
set sel [atomselect top "segname $templateSeg"]
set res0List [lsort -integer -unique -index 0 [$sel get {resid resname}]]
$sel delete

set resList {}
foreach resItem $res0List {
    set res [lindex $resItem 0]
    set resName0 [lindex $resItem 1] 

    set nuc [string index $sequence $nSeq]
    if {[string equal A $nuc]} {
	set loadName ade
	set resName ADE
    } elseif {[string equal C $nuc]} {
	set loadName cyt
	set resName CYT
    } elseif {[string equal G $nuc]} {
	set loadName gua
	set resName GUA
    } else {
	set loadName thy
	set resName THY
    }

    if {[string equal $resName $resName0]} {
	puts "No change: ${resName0} ${res}"

	# Use the old structure and coordinates.
	set s [atomselect top "segname $templateSeg and resid $res"]
	set posList [$s get {x y z}]
	set resNameList [$s get resname]
	set nameList [$s get name]
	$s delete
	
	# Write the pdb lines.
	foreach pos $posList resName $resNameList name $nameList {
	    puts $out [makePdbLineFull $n $templateSeg $res $name $resName $pos]
	    incr n
	}
	
	lappend resList $resItem
    } else {
	# Mutate the structure.
	puts "Mutating: ${resName0}->${resName} ${res}"
	set tempRes $res
	set tempRot [getBaseBasis $templateSeg $tempRes]
	set tempCen [getBasePos $templateSeg $tempRes]
	
	# Write the pdb lines.
	foreach pos $loadPos($loadName) rec $loadRec($loadName) {
	    set nucCen $loadCen($loadName)
	    set nucRot $loadRot($loadName)
	    
	    set r [vecTransform [matTranspose $nucRot] [vecsub $pos $nucCen]]
	    set r [vecadd [vecTransform $tempRot $r] $tempCen]
	    
	    puts $out [makePdbLine $rec $n $templateSeg $res $r]
	    incr n
	}
	lappend resList [list $res [lindex $loadResName($loadName) 0]]
    }

    incr nSeq
}
close $out
mol delete top

# Have psfgen build the structure.
package require psfgen
topology $topFile 
resetpsf
#psfalias

segment $outSeg {
    first 5TER
    last 3TER
    pdb $outSeg.pdb
}

#puts $resList
foreach record $resList {
    foreach {resId resName} $record { break }
    if {[string equal $resName "THY"] || [string equal $resName "CYT"] } {
	patch DEO1 $outSeg:$resId  
    }
    if {[string equal $resName "ADE"] || [string equal $resName "GUA"] } {
	patch DEO2 $outSeg:$resId  
    }
}
coordpdb ${outSeg}.pdb   
guesscoord

writepsf ${outSeg}1.psf
writepdb ${outSeg}1.pdb

# Map the coordinates to what we did earlier.
# That is, undo psfgen's attempts at matching coordinates.
set tempMol [mol load pdb $outSeg.pdb]
set mapMol [mol load psf ${outSeg}1.psf pdb ${outSeg}1.pdb]

foreach record $resList {
    foreach {resId resName} $record { break }
    mapCoordinates $tempMol $mapMol $outSeg $resId all
}

set all [atomselect $mapMol all]
$all writepsf $psf
$all writepdb $pdb
$all delete
mol delete $tempMol 
mol delete $mapMol 
exit
