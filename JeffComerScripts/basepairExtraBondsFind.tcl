# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) and position (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set segA ADNA
set segB BDNA
set spring 1
# Input:
set psf pore6_b-dna_water91924.psf
set pdb fit_pore6_b-dna1.pdb
# Output:
set outPrefix pairs_b-dna

# Define the hydrogen-bond-forming atoms.
set hydrogenBondLight {{ADE N1 THY H3 1.81} {ADE H61 THY O4 1.87} {GUA H21 CYT O2 1.80} {GUA H1 CYT N3 1.84} {GUA O6 CYT H41 1.80}}
set hydrogenBond {{ADE N1 THY N3 2.90} {ADE N6 THY O4 2.88} {GUA N2 CYT O2 2.81} {GUA N1 CYT N3 2.87} {GUA O6 CYT N4 2.81}}

# Make the hydrogen bond list symmetric by adding reverses.
foreach b $hydrogenBond {
    lappend hydrogenBond [list [lindex $b 2] [lindex $b 3] [lindex $b 0] [lindex $b 1] [lindex $b 4]]
}
# Make list of the first resname.
set hBondResName {}
foreach b $hydrogenBond {
    lappend hBondResName [lindex $b 0]
}

# Open the output files.
set out [open ${outPrefix}${spring}.txt w]

# Load the system.
mol load psf $psf pdb $pdb

# Find the complementary bases.
set selA [atomselect top "segname $segA and name C1'"]
set selB [atomselect top "segname $segB and name C1'"]
set resAList [lsort -index 0 -real [$selA get {z resid}]]
set resBList [lsort -index 0 -real [$selB get {z resid}]]

set resA {}
set resB {}
foreach a $resAList b $resBList {
    lappend resA [lindex $a 1]
    lappend resB [lindex $b 1]
}

# Get the residue names.
set resNameA {}
set resNameB {}
foreach a $resA b $resB {
    set selA [atomselect top "segname $segA and resid $a"]
    set selB [atomselect top "segname $segB and resid $b"]
    lappend resNameA [lindex [$selA get resname] 0]
    lappend resNameB [lindex [$selB get resname] 0]
    $selA delete
    $selB delete
}

# Determine the indices of the atoms forming the base hydrogen bonds.
set hbA {}
set hbB {}
foreach a $resA b $resB rnA $resNameA rnB $resNameB {
    set indList [lsearch -all $hBondResName $rnA]
    
    # Find the H-bond forming atoms for these two residues.
    set nBonds 0
    foreach ind $indList {
	foreach {baseA nameA baseB nameB length0} [lindex $hydrogenBond $ind] {break}
	
	if {[string equal $baseB $rnB]} {
	    set selA [atomselect top "segname $segA and resid $a and name $nameA"]
	    set selB [atomselect top "segname $segB and resid $b and name $nameB"]
	    
	    if {[$selA num] != 1 && [$selB num] != 1}  {
		puts stderr "Warning! Cannot bond $baseA to $baseB."
		puts stderr "Atoms $a:$nameA and $b:$nameB were not found."
	    } else {
		set iA [lindex [$selA get index] 0]
		set iB [lindex [$selB get index] 0]
		lappend hbA $iA
		lappend hbB $iB
		incr nBonds
		puts "bond: $a $rnA $b $rnB"
		puts $out "bond $iA $iB $spring $length0"
	    }
	    $selA delete
	    $selB delete
	}
    }

    if {$nBonds == 0} {
	puts stderr "Warning! Found no bonds between residues $a $rnA and $b $rnB."
    } else {
	puts $out ""
    }
}
puts "Found [llength $hbA] potential hydrogen bonds!\n"


close $out
exit



