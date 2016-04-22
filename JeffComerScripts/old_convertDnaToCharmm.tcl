# vmd -dispdev text -e convertDNAtoCHARMM.tcl
set topFile top_all27_prot_na.inp
set segName BDNA
# Input:
set psf0 longDNA_AMBER.psf
set pdb0 longDNA_AMBER.pdb
# Output:
set psf longDNA_ADNA.psf
set pdb longDNA_ADNA.pdb

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

package require psfgen
topology $topFile
resetpsf
psfalias

########################################################################
readpsf  $psf0
coordpdb $pdb0

mol load psf $psf0 pdb $pdb0

set selDNA [atomselect top "segid $segName"] 
$selDNA writepdb $segName.pdb
set tmpList [lsort -unique [$selDNA get {resid resname}] ]
delatom $segName

writepsf rest.psf
writepdb rest.pdb



resetpsf


readpsf  rest.psf
coordpdb rest.pdb 


segment $segName {
    first 5TER
    last 3TER
    pdb $segName.pdb
}

puts $tmpList
foreach record $tmpList {
    foreach {resid resname} $record { break }
    if {$resname == "THY" || $resname == "CYT" } {
	patch DEO1 $segName:$resid  
    }
    if {$resname == "ADE" || $resname == "GUA" } {
	patch DEO2 $segName:$resid  
    }
}
coordpdb $segName.pdb   


guesscoord

writepsf $psf
writepdb $pdb

mol delete top
exit

