# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# jcomer2@uiuc.edu

# Parameters:
set simNum 20
set selText "water or ions"
set coorFileList [glob /scratch/tbgl/jcomer/basepair/output/b*_gc_pos*.coor]
# Input:
set inName pore_gcm_basepair
set refName pore_gc_basepair
# Output:
set outPrefix ../init/init_gcm_pos

# Load the system.
set refMol [mol load psf $refName.psf pdb $refName.pdb]
set inMol [mol load psf $inName.psf pdb $inName.pdb]
set coorNum [llength $coorFileList]

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

# Make the selections.
set sel [atomselect $inMol $selText]
set all [atomselect $inMol all]
set refSel [atomselect $refMol $selText]

# Write the smd files.
for {set sim 0} {$sim < $simNum} {incr sim} {
    puts "SIM $sim"

    # Set the coordinates from a random coordinate file.
    set coorInd [expr {int(floor(rand()*$coorNum))}]
    set coorFile [lindex $coorFileList $coorInd] 
    puts "USING [$refSel num] coordinates from `$coorFile'."
    animate delete beg 0 end -1 $refMol
    animate read namdbin $coorFile waitfor all $refMol

    set posList [$refSel get {x y z}]
    puts [lindex $posList 0]
    $sel set {x y z} $posList

    $all writepdb ${outPrefix}${sim}.pdb
}

$sel delete
$all delete
$refSel delete
mol delete top
exit
