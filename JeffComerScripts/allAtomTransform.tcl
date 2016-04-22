# Author: Jeff Comer <jcomer2@illinois.edu>

# Input:
set transformFileList [glob *.trans]
set inSegName "ADNA"
set outSegName "D"
set temp tmp

set gridList {} 
# resName matchString structPrefix
lappend gridList [list ADE *pot* ../../4_canonical/typical_ade-thy]
#lappend gridList [list ADE *ade_* nucleo_ade_pot]
#lappend gridList [list THY *thy_* nucleo_thy_pot]
#lappend gridList [list GUA *gua_* nucleo_gua_pot]
#lappend gridList [list CYT *cyt_* nucleo_cyt_pot]

package require psfgen
proc silent {} {}

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	set t [lindex $tok 0]
	lappend r $tok
    }

    close $in
    return $r
}

proc matMake4 {m {d {0.0 0.0 0.0}}} {
    set ret {}
    lappend ret [concat [lindex $m 0] [lindex $d 0]]
    lappend ret [concat [lindex $m 1] [lindex $d 1]]
    lappend ret [concat [lindex $m 2] [lindex $d 2]]
    lappend ret [list 0.0 0.0 0.0 1.0]
}

# Load the molecules.
set molList {}
set selList {}
set posList {}
foreach g $gridList {
    foreach {resName matchString structPrefix} $g { break }
    lappend molList [mol load psf $structPrefix.psf pdb $structPrefix.pdb]
    set sel [atomselect top "segname $inSegName"]
    lappend selList $sel
    lappend posList [$sel get {x y z}]
}

foreach transformFile $transformFileList {
    set data [readData $transformFile]

    set indexList {}
    set transList {}
    foreach d $data {
	set grid [lindex $d 0]
	set trans {}

	# Find the structure corresponding to the grid.
	set gInd 0
	foreach gridSet $gridList {
	    set gridMatch [lindex $gridSet 1]
	    if {[string match "$gridMatch" $grid]} {
		lappend prefixList [lindex $gridSet 2]
		set rot [list [lrange $d 1 3] [lrange $d 4 6] [lrange $d 7 9]]
		set disp [lrange $d 10 12] 

		#puts "rot: $rot"
		#puts "disp: $disp"

		lappend transList [matMake4 $rot $disp]
		lappend indexList $gInd
		puts "$grid -> [lindex $gridSet 2]"
		break
	    }
	    incr gInd
	}
    }

    # Combine the structures.
    resetpsf

    set count 0
    # Do the transformations.
    foreach index $indexList trans $transList {
	set sel [lindex $selList $index]
	set pos [lindex $posList $index]

	# Reset the positions.
	$sel set {x y z} $pos

	# Set the new segment name.
	$sel set segname ${outSegName}${count}

	# Transform to the new positions.
	puts $trans
	$sel move $trans
	$sel writepsf $temp.psf
	$sel writepdb $temp.pdb

	# Combine with psfgen.
	readpsf $temp.psf
	coordpdb $temp.pdb

	incr count
    }

    # Write the final structure.
    writepsf $transformFile.psf
    writepdb $transformFile.pdb
}

exit
