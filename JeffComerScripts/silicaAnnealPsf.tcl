# Create a silica psf from a pdb with correct type names,
# correct charges, and no bonds.
# Use with: vmd -dispdev text -e silicaAnnealPsf.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set pdbFile silica_shell_ready1.pdb
# Output:
set psfFile silica_shell_ready1.psf
# Si parameters
set nameSi "SI.*"
set massSi 28.085500
set chargeSi 2.4
set typeSi SI
# O parameters
set nameO "O.*"
set massO 15.9994
# The charge will be set for neutrality later.
set chargeO -1.2
set typeO O

proc main {} {
	global pdbFile psfFile
	global nameSi chargeSi nameO chargeO
	
	# Load the pdb.
	mol load pdb $pdbFile
	set nAtoms [molinfo top get numatoms]
	
	# Compute the charge on O to enforce neutrality.
	set sel [atomselect top "name \"${nameSi}\""]
	set nSi [$sel num]
	$sel delete
	set sel [atomselect top "name \"${nameO}\""]
	set nO [$sel num]
	$sel delete
	mol delete all
	set chargeO0 $chargeO
	set chargeO [expr -$chargeSi*$nSi/$nO]
	puts "\nSet charge on O to $chargeO to enforce neutrality."
	puts "Relative error: [expr ($chargeO-$chargeO0)/$chargeO0]\n"
			
	puts "Writing psf file..."
	manifestPsf $psfFile $pdbFile $nAtoms
	puts "The file $psfFile was written successfully."
}

# Write the psf file.
proc manifestPsf {psfFile pdbFile nAtoms} {
	global nameSi massSi chargeSi typeSi
	global nameO massO chargeO typeO
	 
	set dummy "          0"
	set totalBonds 0
	set totalAngles 0
	set out [open $psfFile w]
	
	##### HEADER
	puts $out "PSF"
	puts $out ""
	puts $out "       1 !NTITLE"
	puts $out " REMARKS original generated structure x-plor psf file"

	##### ATOMS
	puts "Writing atom records..."
	puts $out ""
	puts $out "[format %8i $nAtoms] !NATOM"
	
	# Open the pdb to extract the atom records.
	set inStream [open $pdbFile r]
	set atom 1
	foreach line [split [read $inStream] \n] {
    		set string0 [string range $line 0 3]
		if {![string match $string0 "ATOM"]} {continue}
		
		# Extract each pdb field.
		set record [string range $line 0 5]
		set serial [string range $line 6 10]
		set name [string range $line 12 15]
		set altLoc [string range $line 16 16]
		set resName [string range $line 17 19]
		set chainId [string range $line 21 21]
		set resId [string range $line 22 25]
		set iCode [string range $line 26 26]
		set x [string range $line 30 37]
		set y [string range $line 38 45]
		set z [string range $line 46 53]
		set occupancy [string range $line 54 59]
		set beta [string range $line 60 65]
		set segName [string range $line 72 75]
		set element [string range $line 76 77]
		set charge [string range $line 78 79]
		
		# Write the atom record.	
		puts -nonewline $out [format "%8i " $atom]
		puts -nonewline $out [format "%-4s " $segName]
		puts -nonewline $out [format "%-4i " $resId]
		puts -nonewline $out [format "%-3s " $resName]
		puts -nonewline $out [format "%-4s  " $name]
		if {[regexp $nameSi $name]} {
			puts -nonewline $out [format "%-4s  " $typeSi]
			puts -nonewline $out [format "% 5.6f       " $chargeSi]
			puts -nonewline $out [format "%6.4f " $massSi]
		} else {
			puts -nonewline $out [format "%-4s  " $typeO]
			puts -nonewline $out [format "% 5.6f       " $chargeO]
			puts -nonewline $out [format "%6.4f " $massO]
		}
		puts $out $dummy
		
		incr atom
	}
	close $inStream
	puts $out ""
  
	##### BONDS
	# Write the bonds.
	set total [format %8i $totalBonds]
	puts $out "$total !NBOND: bonds"
	set num 0
	puts $out ""

	##### ANGLES
	# Write the angles.
	puts $out "[format %8i $totalAngles] !NTHETA: angles"
	set num 0
	puts $out ""

	# Write everything else.
	##### DIHEDRALS
	set nDihedrals 0
	puts $out ""
	puts $out "[format %8i $nDihedrals] !NPHI: dihedrals"
	puts $out ""

	##### IMPROPERS
	set nImpropers 0
	puts $out ""
	puts $out "[format %8i $nImpropers] !NIMPHI: impropers"
	puts $out ""

	##### DONORS
	set nDonors 0
	puts $out ""
	puts $out "[format %8i $nDonors] !NDON: donors"
	puts $out ""

	##### ACCEPTORS
	set nAcceptors 0
	puts $out ""
	puts $out "[format %8i $nAcceptors] !NACC: acceptors"
	puts $out ""

	##### NON-BONDED
	set nNB 0
	puts $out ""
	puts $out "[format %8i $nNB] !NNB"
	puts $out ""

	set tmp [expr int($nAtoms/8)]
	set tmp2 [expr $nAtoms -$tmp*8]
	for {set i 0} {$i <$tmp} {incr i} {
		puts $out "       0       0       0       0       0       0       0       0"
	}
	set lastString ""
	for {set i 0} {$i <$tmp2} {incr i} {
	    set lastString "${lastString}       0"
	}
	puts $out $lastString

	####### GROUPS
	puts $out ""
	puts $out "       1       0 !NGRP"
	puts $out "       0       0       0"
	puts $out ""
	puts $out ""
	close $out
}

main
exit



