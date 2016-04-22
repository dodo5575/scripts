# Use with: vmd -dispdev text -e convertDnaToCharmm.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set nameList {DNA.170 DNA.320 DNA.1000 CG.170}
set topFile top_all27_na.rtf
set segList {ADNA BDNA}
# Input:
# Output:
set otherName other

proc psfalias {}  {
	# Define common aliases
	# Here's for nucleics
	alias residue G GUA
	alias residue C CYT
	alias residue A ADE
	alias residue T THY
	alias residue U URA

	foreach bp { GUA CYT ADE THY URA } {
		alias atom $bp "O5\*" O5'
		alias atom $bp "C5\*" C5'
		alias atom $bp "O4\*" O4'
		alias atom $bp "C4\*" C4'
		alias atom $bp "C3\*" C3'
		alias atom $bp "O3\*" O3'
		alias atom $bp "C2\*" C2'
		alias atom $bp "O2\*" O2'
		alias atom $bp "C1\*" C1'
	}

	alias atom ILE CD1 CD
	alias atom SER HG HG1
	alias residue HIS HSE

	# Heme aliases
	alias residue HEM HEME
	alias atom HEME "N A" NA
	alias atom HEME "N B" NB
	alias atom HEME "N C" NC
	alias atom HEME "N D" ND

	# Water aliases
	alias residue HOH TIP3
	alias atom TIP3 O OH2

#	alias atom HOH O OH2

	# Ion aliases
	alias residue NA SOD
	alias atom NA NA SOD
	alias residue K POT
	alias atom K K POT
	alias residue ICL CLA
	alias atom ICL CL CLA
}

########################################################################
package require psfgen
topology $topFile

foreach name $nameList {
    set psf ${name}.psf
    set pdb ${name}.pdb
    set outName charmm_${name}

    resetpsf
    psfalias

    mol load psf $psf pdb $pdb
    set all [atomselect top all]
    set allN [$all num]
    $all delete

    # Change the atom names.
    set nameChange {{"resname THY and name C7" C5M} {"resname THY and name H71" H51} {"resname THY and name H72" H52} {"resname THY and name H73" H53}}
    foreach change $nameChange {
	foreach {text name} $change {break}
	set s [atomselect top $text]
	$s set name $name
	puts "Changed the name of [$s num] atoms defined by ($text) to $name."
	$s delete
    }

    # Read the structure into psfgen.
    readpsf $psf
    coordpdb $pdb

    # Write the DNA segments.
    set segResList {}
    set n 0
    foreach seg $segList {
	set sel [atomselect top "segname $seg"]
	$sel writepdb $seg.pdb
	set n [expr $n + [$sel num]]
	lappend segResList [lsort -unique -index 0 -integer [$sel get {resid resname}]]
	$sel delete

	# Delete the DNA from psfgen.
	delatom $seg
    }

    # Write the everything else that's not DNA.
    writepsf $otherName.psf
    writepdb $otherName.pdb

    # Build the DNA.
    resetpsf
    foreach seg $segList resList $segResList {
	segment $seg {
	    first none
	    last none
	    pdb $seg.pdb
	}
	
	foreach res $resList {
	    foreach {resid resname} $res {break}
	    if {[string equal $resname "THY"] || [string equal $resname "CYT"]} {
		patch DEO1 $seg:$resid  
	    }
	    if {[string equal $resname "ADE"] || [string equal $resname "GUA"]} {
		patch DEO2 $seg:$resid  
	    }
	}
	coordpdb $seg.pdb

	# Link to make effectively infinite DNA.
	set resFirst [lindex $resList 0 0]
	set resLast [lindex $resList end 0]
	puts "Linking $seg:$resFirst $seg:$resLast"
	patch LKNA $seg:$resLast $seg:$resFirst 
	regenerate angles dihedrals
    }

    # Write the completed DNA structures.
    guesscoord
    writepsf ${outName}.0.psf
    writepdb ${outName}.0.pdb

    # Combine the DNA with everything else.
    resetpsf
    readpsf ${outName}.0.psf
    coordpdb ${outName}.0.pdb
    if {$n < $allN} {
	readpsf $otherName.psf
	coordpdb $otherName.pdb
	writepsf ${outName}.psf
	writepdb ${outName}.pdb
    }

    mol delete top
}
exit
