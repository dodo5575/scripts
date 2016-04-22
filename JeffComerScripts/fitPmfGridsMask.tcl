# Match one selection with another.
# Use with: vmd -dispdev text -e fit.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

#set bp at
#set ion pot
set selText "segname ADNA BDNA"
# Input:
set baseGrid ener92_${bp}_${ion}.dx
set refPsf pore_${bp}_basepair.psf
set refCoor pore_${bp}_basepair.pdb
# Contains "resname" "grid file" "corresponding structure prefix" for each type of residue
set fitText "segname ADNA"
set gridList {} 
lappend gridList [list ADE pmf_pmfHard_ade_${ion}.dx mask_ade_${ion}.dx nucleo_ade_pot]
lappend gridList [list THY pmf_pmfHard_thy_${ion}.dx mask_thy_${ion}.dx nucleo_thy_pot]
lappend gridList [list GUA pmf_pmfHard_gua_${ion}.dx mask_gua_${ion}.dx nucleo_gua_pot]
lappend gridList [list CYT pmf_pmfHard_cyt_${ion}.dx mask_cyt_${ion}.dx nucleo_cyt_pot]
# Outfile:
set outFile trans_${bp}_${ion}.txt

proc getResFit {fitMol fitText refMol refText} {
    set fitSel [atomselect $fitMol $fitText]
    set fitList [$fitSel get {name resname}]
    set refSel [atomselect $refMol $refText]
    set refList [$refSel get {name resname}]

    # Filter atoms that aren't in both selections.
    set masterList {}
    foreach resSet $fitList {
	set ind [lsearch $refList $resSet]
	if {$ind >= 0}  {
	    lappend masterList $resSet
	} else {
	    #puts "No $resSet."
	}
    }
    if {[llength $masterList] < 3} {
	puts "Not enough atoms ([llength $masterList]) for comparison! `$fitText' `$refText'"
	exit
    }
       
    set masterText "(name [lindex $masterList 0 0] and resname [lindex $masterList 0 1])"
    for {set i 1} {$i < [llength $masterList]} {incr i} {
	set masterText "$masterText or (name [lindex $masterList $i 0] and resname [lindex $masterList $i 1])"
    }
    set fitSel1 [atomselect $fitMol "($fitText) and ($masterText)"]
    set refSel1 [atomselect $refMol "($refText) and ($masterText)"]
    set fitList1 [lsort [$fitSel1 get {name resname}]]
    set refList1 [lsort [$refSel1 get {name resname}]]

    # Get the order list.
    set order {}
    foreach name $fitList1 {
	set ind [lsearch $refList1 $name]
	if {$ind < 0 } {
	    puts $refList
	    puts "Cannot find atom $name in reference molecule ${refText}!"
	    exit
	}
	lappend order $ind
    }
    
    set ret [measure fit $fitSel1 $refSel1 weight mass order $order]
    $fitSel delete
    $refSel delete
    $fitSel1 delete
    $refSel1 delete
    return $ret
}

proc oneMatrixToAnother {m} {
    set row0 [lrange [lindex $m 0] 0 2]
    set row1 [lrange [lindex $m 1] 0 2]
    set row2 [lrange [lindex $m 2] 0 2]
    set disp [list [lindex $m 0 3] [lindex $m 1 3] [lindex $m 2 3]]
    return [list $row0 $row1 $row2 $disp]
}

# Open the output.
set out [open $outFile w]
puts $out "$baseGrid $baseGrid 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0"

# Load the molecules.
set refMol [mol load psf $refPsf]
mol addfile $refCoor

set molList {}
foreach g $gridList {
    foreach {resName gridFile maskFile structPrefix} $g { break }
    lappend molList [mol load psf $structPrefix.psf pdb $structPrefix.pdb]
}

# Get the residues.
set sel [atomselect $refMol $selText]
set resList [lsort -unique [$sel get {segname resid resname}]]
$sel delete

# Loop through the residues of the selection.
foreach r $resList {
    # Get the resName.
    foreach {seg res resName} $r { break }
    set refText "segname $seg and resid $res"

    # Find the corresponding grid.
    set chosen -1
    set j 0
    foreach g $gridList {
	if {[string equal [lindex $g 0] $resName]} {
	    set chosen $j
	    break
	}
	incr j
    }

    # Check that we found the grid.
    if {$chosen < 0} {
	puts "WARNING! No grid was found for resname `$resName'."
	continue
    }
    puts "RESIDUE: $resName"
    
    # Get the transformation matrix.
    set gridFile [lindex $gridList $chosen 1]
    set maskFile [lindex $gridList $chosen 2]
    set fitMol [lindex $molList $chosen]
    #puts "REF: $fitMol $fitText $refMol $refText"
    set basis [getResFit $fitMol $fitText $refMol $refText]
    
    # Write the results.
    set outText "$gridFile $maskFile [join [oneMatrixToAnother $basis]]"
    puts $out $outText
    #puts $outText
}
close $out

mol delete $refMol
foreach m $molList {
    mol delete $m
}

#exit
