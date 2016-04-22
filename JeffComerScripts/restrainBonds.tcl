# Add harmonic constraints to silicon nitride.
# to use: vmd -dispdev text -e constrainSilicon.tcl
# jcomer2@uiuc.edu

# Parameters:
set fixedText "index 0"
set restText "index 1"
set name chl-chl
set psf neu_chl-chl.psf
set pdb neu_chl-chl.pdb
set occupancy 1.0
set stride 1
# Input:
set indexFile win_ions.txt
# Output:
set outDir bonds

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

mol load psf $psf pdb $pdb
set all [atomselect top all]
set fixedSel [atomselect top $fixedText]
set restSel [atomselect top $restText]
$fixedSel set {x y z} {{0 0 0}}

# Read the index file.
set data [readData $indexFile]
set nodePos {}
set nodeSpring {}
set nodeIndex {}
foreach d $data {
    lappend nodePos [lrange $d 1 3]
    lappend nodeSpring [lrange $d 4 6]
    lappend nodeIndex [lindex $d 0]
}

foreach pos $nodePos spring $nodeSpring index $nodeIndex {
    set outFile $outDir/${name}_pos${index}.pdb

    set k [expr 2.0*[lindex $spring 2]]
    set d [lindex $pos 2]
    
    $restSel set {x y z} [list [list 0 0 $d]]
    $all set beta 0.0
    $restSel set beta $k

    $all writepdb $outDir/restrain_${name}${index}.pdb
}

$restSel delete
$fixedSel delete
$all delete
mol delete all
exit
