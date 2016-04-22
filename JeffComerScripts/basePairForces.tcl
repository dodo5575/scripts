# This script makes a file for the H-bonds between nitrogen bases.
# The file is used by the TCL forces script confinePore+Pairs.tcl.
# Use with: vmd -dispdev text -e restrainBasePairs.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set dnapdb "dsDNA-all.pdb"
set dnapsf "dsDNA-all.psf"
set dnabincoord "equil-DS.restart.coor"
# Output:
set namdConstFile "basePairs.txt"

# Load the DNA.
mol load psf $dnapsf pdb $dnapdb
mol addfile $dnabincoord waitfor all

# Open the constraint file and write the header.
set out [open $namdConstFile w]

# Set the spring constants of the restraints.
set kf1 [expr 10.0*0.2]
set kf2 [expr 20.0*0.2]
set kf3 [expr 15.0*0.2]

# Loop through the nucleotides of strand A.
set selA [atomselect top "segid ADNA"]
set reslist [lsort -unique -integer [$selA get resid]]
foreach resA $reslist {
	# Get the name of the nucleotide.
	set resSelA [atomselect top "segid ADNA and resid $resA"]
	set nucleotideA [lindex [$resSelA get resname] 0]

	# Find the position of an atom on the nucleotide that interfaces with the other strand.
	if {$nucleotideA == "ADE"} {
		set siteA [atomselect top "segid ADNA and resid $resA and name N1"]
	} elseif {$nucleotideA == "THY"} {
		set siteA [atomselect top "segid ADNA and resid $resA and name N3"]
	} elseif {$nucleotideA == "GUA"} {
		set siteA [atomselect top "segid ADNA and resid $resA and name N1"]
	} elseif {$nucleotideA == "CYT"} {
		set siteA [atomselect top "segid ADNA and resid $resA and name N3"]
	} else {
		puts "ERROR: Unknown nucleotide $nucleotideA"
		close $out
		exit
	}
	set ra [lindex [$siteA get {x y z}] 0]

	# First find atoms on strand B within 5.0A of this nucleotide.
	set selB [atomselect top \
	"segid BDNA and within 5.0 of (segid ADNA and resid $resA)"]
	
	# Identify the atom on strand B closest to the appropriate atom on strand A.
	set rb [$selB get {resid x y z}]
	set db {}
	foreach r $rb {
		set dist2 [veclength2 [vecsub [lrange $r 1 3] $ra]]
		lappend db [list [lindex $r 0] $dist2]
	}
	set db [lsort -real -index 1 $db]
	
	# The matching nucleotide on strand B is the residue of the closest atom.
	set resB [lindex $db 0 0]
	set resSelB [atomselect top "segid BDNA and resid $resB"]
	
	
	set nucleotideB [lindex [$resSelB get resname] 0]
	puts "Pair: ($resA, $nucleotideA) ($resB, $nucleotideB)"
	
	# Write the constraints in the format "segname resid atomname segname resid atomname kf ref".
	if {$nucleotideA == "ADE" && $nucleotideB == "THY" } {
		puts $out "ADNA $resA N6 BDNA $resB O4 $kf1 2.9"
		puts $out "ADNA $resA N1 BDNA $resB N3 $kf2 2.88"
		#puts $out "ADNA $resA H2 BDNA $resB O2 $kf3 2.95"
	} elseif {$nucleotideA == "THY" && $nucleotideB == "ADE" } {
		puts $out "ADNA $resA O4 BDNA $resB N6 $kf1 2.9"
		puts $out "ADNA $resA N3 BDNA $resB N1 $kf2 2.88"
		#puts $out "ADNA $resA O2 BDNA $resB H2 $kf3 2.95"
	} elseif {$nucleotideA == "GUA" && $nucleotideB == "CYT" } {
		puts $out "ADNA $resA O6 BDNA $resB N4 $kf1 2.81"
		puts $out "ADNA $resA N1 BDNA $resB N3 $kf2 2.87"
		puts $out "ADNA $resA N2 BDNA $resB O2 $kf3 2.81"
	} elseif {$nucleotideA == "CYT" && $nucleotideB == "GUA" } {
		puts $out "ADNA $resA N4 BDNA $resB O6 $kf1 2.81"
		puts $out "ADNA $resA N3 BDNA $resB N1 $kf2 2.87"
		puts $out "ADNA $resA O2 BDNA $resB N2 $kf3 2.81"
	} else {
		puts "ERROR: Unknown pairing of nucleotides"
		close $out
		exit
	}
}
close $out  
exit




