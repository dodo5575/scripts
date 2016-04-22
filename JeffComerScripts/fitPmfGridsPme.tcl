# Match one selection with another.
# Use with: vmd -dispdev text -e fit.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set ionList {pot}
set selText "nucleic"
# Input:
set pdbList [glob ../struct/*_pot.pdb]
set fitText "segname ADNA"
# Outfile:
set outDir transform

proc main {pdb ion} {
    global selText fitText outDir
    set name [trimLastUnder [trimFirstUnder [trimPath [trimExtension $pdb]]]]

    set baseGrid ../struct/helix.dx

    # Make a list of corresponding grids, point charge files, and all-atoms structures.
    # Contains "resname" "grid file" "point charge file" "corresponding structure prefix" for each type of residue
    set gridList {} 
    lappend gridList [list ADE prepare/pmf_readyLow_ade_${ion}.dx nucleo/nucleo_ade_pot]
    lappend gridList [list THY prepare/pmf_readyLow_thy_${ion}.dx nucleo/nucleo_thy_pot]
    #lappend gridList [list GUA prepare/pmf_readyLow_gua_${ion}.dx nucleo/nucleo_gua_pot]
    #lappend gridList [list CYT prepare/pmf_readyLow_cyt_${ion}.dx  nucleo/nucleo_cyt_pot]
    #lappend gridList [list MCY prepare/pmf_readyLow_mcyt_${ion}.dx nucleo/nucleo_mcyt_pot]

    if {![file exists $baseGrid]} {
	puts "$baseGrid does not exist."
	exit
    }

    set refPdb $pdb 
    set refName [trimLastUnder [trimExtension [trimPath $refPdb]]]
    set outFile $outDir/${refName}_${ion}.trans

    # Open the output.
    set out [open $outFile w]
    puts $out "$baseGrid 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0"

    # Load the molecules.
    set refMol [mol load pdb $refPdb]

    set molList {}
    foreach g $gridList {
	foreach {resName gridFile structPrefix} $g { break }
	lappend molList [mol load psf $structPrefix.psf pdb $structPrefix.pdb]
    }

    # Get the residues.
    set sel [atomselect $refMol $selText]
    set resList [lsort -unique -index 0 -integer [$sel get {residue resname}]]
    $sel delete

    # Loop through the residues of the selection.
    foreach r $resList {
	# Get the resName.
	foreach {residue resName} $r { break }
	set refText "residue $residue"

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
	#set chargeFile [lindex $gridList $chosen 2]
	set fitMol [lindex $molList $chosen]
	#puts "REF: $fitMol $fitText $refMol $refText"
	set basis [getResFit $fitMol $fitText $refMol $refText]
	
	# Write the results.
	#set outText "$gridFile $chargeFile [join [oneMatrixToAnother $basis]]"
	set outText "$gridFile [join [oneMatrixToAnother $basis]]"
	puts $out $outText
	#puts $outText
    }
    close $out

    mol delete $refMol
    foreach m $molList {
	mol delete $m
    }
}

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
    set fitList1 [$fitSel1 get {name resname}]
    set refList1 [$refSel1 get {name resname}]

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

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

proc trimFirstUnder {name} {
    set ind [string first "_" $name]
    return [string range $name [expr $ind+1] end]
}

proc trimLastUnder {name} {
    set ind [string last "_" $name]
    return [string range $name 0 [expr $ind-1]]
}

foreach ion $ionList {
    foreach pdb $pdbList {
	    main $pdb $ion
    }
}

exit
