# Use with: vmd -dispdev text -e generateHairpin.tcl
# Author: Jeff Comer <jcomer2@uiuc.edu>

set range 4
set substrateText "resname SIO2"
set selText "type SI and within $range of name OH2"
set collideName N1
set collideRange 1.5
set outSeg AMI
# Input:
set topFile top_amine_silicon_charge.inp
set topFile1 top_all27_prot_lipid_pot.inp
set psf trap2.2_KCl.psf
set pdb silica_trap2.2_0Va3.pdb
set templateStruct one_amine
set densityGrid pore_density3.dx
# Output:
set outName amine_coating
set substrateName substrate

source vector.tcl
set templatePsf $templateStruct.psf
set templatePdb $templateStruct.pdb

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

# Define a basis for the nucleotide.
proc getAmineBasis {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
    
    # Get the amine basis.
    set selX0 [atomselect $mole "($selText) and name SI1"]
    set selX1 [atomselect $mole "($selText) and name N1"]
    set selY0 [atomselect $mole "($selText) and name SI1"]
    set selY1 [atomselect $mole "($selText) and name C4"]

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

# Define a basis for the nucleotide.
proc getAminePos {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
    
    # Get the amine basis.
    set sel [atomselect $mole "($selText) and name SI1"]
    return [lindex [$sel get {x y z}] 0]
}

proc deleteAtoms {atomList inPsf inPdb outName} {
    resetpsf
    readpsf $inPsf
    coordpdb $inPdb

    # Delete the selection.
    foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
    }

    writepsf $outName.psf
    writepdb $outName.pdb
    puts "[llength $atomList] atoms were deleted."
    return
}

proc combine {inName0 inName1 outName} {
    resetpsf
    
    readpsf $inName0.psf
    coordpdb $inName0.pdb
    readpsf $inName1.psf
    coordpdb $inName1.pdb

    writepsf $outName.psf
    writepdb $outName.pdb

    return
}

set templateRec [extractPdbRecords $templatePdb]
set templateCoord [extractPdbCoords $templatePdb]

# Load the template and get its info.
mol load psf $templatePsf pdb $templatePdb
set sel [atomselect top all]
set addSeg [lindex [$sel get segname] 0]
set addRes [lindex [$sel get resid] 0]
set tempBasis [getAmineBasis $addSeg $addRes]
set tempPos [getAminePos $addSeg $addRes]
$sel delete
mol delete top

mol load psf $psf pdb $pdb
set sel [atomselect top $selText]
set subNum [$sel num]
foreach zero {0} {set subPos [$sel get {x y z}]}
$sel delete
set other [atomselect top "(not ($substrateText)) or ($selText)"]
foreach zero {0} {set otherAtoms [$other get {segname resid name}]}
$other delete

puts "Building the PDB."
puts "Chemically modifying $subNum sites."
set out [open $outSeg.pdb w]
set index 1
set res 1
foreach p $subPos {
    foreach pos $templateCoord rec $templateRec {
	set r [vecadd $p [vecsub $pos $tempPos]]
	set line [makePdbLine $rec $index $outSeg $res $r]
	puts $out $line

	incr index
    }
    incr res
}
close $out

# Have psfgen build the structure.
puts "Building the PSF with psfgen."
package require psfgen
topology $topFile 
topology $topFile1
resetpsf
#psfalias

segment $outSeg {
    first NONE
    last NONE
    pdb $outSeg.pdb
}

coordpdb $outSeg.pdb
guesscoord

writepsf ${outName}0.psf
writepdb ${outName}0.pdb

# Delete the coated atoms.
puts "Deleting the coated atoms from the original system."
deleteAtoms $otherAtoms $psf $pdb $substrateName

# Combine the resultant systems.
puts "Combining the systems."
combine $substrateName ${outName}0 ${outName}1

# Rotate the substituents to avoid clashes.
puts "Moving the substituents to avoid clashes."
mol load psf ${outName}1.psf pdb ${outName}1.pdb
set res 1
set maxTries 40
foreach p $subPos {
    set s [atomselect top "segname $outSeg and resid $res and name $collideName and within $collideRange of not (segname $outSeg and resid $res)"]
    set try 0

    while {[$s num] > 0} {
	$s delete

	set s1 [atomselect top "segname $outSeg and resid $res"]
	set rot [matMake4 [matRandomRot]]
	$s1 moveby [vecinvert $p]
	$s1 move $rot
	$s1 moveby $p
	$s1 delete
 
	set s [atomselect top "segname $outSeg and resid $res and name $collideName and within $collideRange of not (segname $outSeg and resid $res)"]
	incr try
	if {$try >= $maxTries} {
	    puts "Could not resolve clash, residue $res"
	    break
	}
    }
    incr res
}

set all [atomselect top all]
$all writepdb $outName.pdb
$all writepsf $outName.psf

exit
