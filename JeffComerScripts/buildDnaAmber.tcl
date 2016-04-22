# this script take a template pdf/pdb files and mutate basis of the nucleic
# acids according to the specified sequence.
# This script makes double stranded DNA only.
# Created by Aleksei Aksimentiev alek@ks.uiuc.edu
# usage:  vmd -dispdev text -e buildDnaAmber.tcl

# Parameters:
# THIS IS THE DNA SEQUENCE
set sequenceDNA ""
# THIS IS THE DNA LENGTH
set sequenceN 24
# to mutate dna according to the $sequence put "yes" 
# otherwise put "no" and specify $sequenceN
set mutate no

# Input:
set templatePdb longDNA_AMBER.pdb
set templatePsf longDNA_AMBER.psf
set topoFile cornell.rtf
# Output:
set newfile dsDnaAmber

####################################################################
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
topology $topoFile
resetpsf
psfalias

########################################################################

## Begining of the main part

set newpdb  [file join $newfile ${newfile}.pdb]
set newpsf  [file join $newfile ${newfile}.psf]
set backbonepdb [file join $newfile ${newfile}-backbone.pdb]
set namdConst [file join $newfile ${newfile}-const.namd]

file mkdir $newfile

readpsf  $templatePsf
coordpdb $templatePdb


set inputID [mol load psf $templatePsf pdb $templatePdb] 
set sel [atomselect top all]
set maxResid [llength [lsort -unique [$sel get resid]]]
puts "makResid $maxResid"
if { $sequenceN > $maxResid} {
    puts "ERROR: the specified DNA sequence is longer than the template"
    puts "EXIT"
    exit
}

set sel [atomselect $inputID "segid ADNA and resid 1 to $sequenceN"]
set listA [lsort -unique [$sel get {resid resname}] ]
$sel writepdb ADNA.pdb


set sel [atomselect $inputID "segid BDNA and resid [expr $maxResid -$sequenceN +1] to $maxResid"]
set listB [lsort  -unique [$sel get {resid resname}] ]
$sel writepdb BDNA.pdb
# splitting the input pdb

set selFirst [atomselect $inputID "(segid ADNA and resid 1 and name O5') or (segid BDNA and resid [expr $maxResid -$sequenceN+1] and name O5')"]
set firstResList [$selFirst get resname]
set firstResidList  [$selFirst get resid]

set selLast [atomselect $inputID "(segid ADNA and resid $sequenceN and name O5') or (segid BDNA and resid $maxResid and name O5')"]     
set lastResList [$selLast get resname]
set lastResidList  [$selLast get resid]

set segidList {ADNA BDNA} 
set patchList [list $listA $listB]

resetpsf

foreach segid $segidList tmpList $patchList firstRes $firstResList lastRes $lastResList firstResid $firstResidList lastResid $lastResidList {

    segment $segid {
	puts "${firstRes}5 ${lastRes}3"
	first ${firstRes}5
	last  ${lastRes}3
#	first 5TER
#	last 3TER
	pdb ${segid}.pdb
    }

    set sel [atomselect top "segid $segid"]
#    set tmpList [$sel get {resid resname}]
#    set tmpList [lsort -unique $tmpList]

    puts $tmpList

    # patch to make DEOXYribose 
    foreach record $tmpList {
	foreach {resid resname} $record { break }

	if {$resname == "ADE"} {
	    if {$resid == $firstResid} {
		patch DOA5 ${segid}:$resid  
	    } elseif {$resid == $lastResid} {
		patch DOA3 ${segid}:$resid  
	    } else {
		patch DOA ${segid}:$resid  
	    }
	}
	if {$resname == "CYT"} {
	    if {$resid == $firstResid} {
		patch DOC5 ${segid}:$resid  
	    } elseif {$resid == $lastResid} {
		patch DOC3 ${segid}:$resid  
	    } else {
		patch DOC ${segid}:$resid  
	    }
	}
	if {$resname == "GUA"} {
	    if {$resid == $firstResid } {
		patch DOG5 ${segid}:$resid  
	    } elseif {$resid == $lastResid } {
		patch DOG3 ${segid}:$resid  
	    } else {
		patch DOG ${segid}:$resid  
	    }
	}
	
	puts "$segid $resid $resname"

    }
    coordpdb ${segid}.pdb   
}

guesscoord

puts "writing built structure to $newpsf $newpdb"

writepsf $newpsf
writepdb $newpdb


mol load psf $newpsf pdb $newpdb 
set sel [atomselect top "resname ADE THY CYT GUA and name C4' O4' C1' C2' P O1P O2P O5' C5' C3' O3' " ]
$sel set beta 1.00
set selALL [atomselect top all]
$selALL writepdb $backbonepdb	 

foreach segid $segidList {
    file delete ${segid}.pdb
}

exit



