 # This script will remove all residues in the selection
# from psf and pdf files.
# Use with: vmd -dispdev text -e removeResidues.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

foreach {dx dy} {18 18 20 20 23 23 24 18 22 17} {
    set atomDiam 0.35; # in nanometers
    set diamX [expr {0.1*$dx + $atomDiam}]; # in nanometers
    set diamY [expr {0.1*$dy + $atomDiam}]; # in nanometers

    set angle 30.0
    set parabLen 1.0; # in nanometers

    set buffer 0.4; # in nanometers
    set frac 0.4; # fraction to delete
    set bufferText "(type NSI and numbonds < 3) or (type SI and numbonds < 4)"

    #Input:
    set inName rough_d9z8_types
    #Output:
    set outName pore${dx}-${dy}

    set psf $inName.psf
    set pdb $inName.pdb
    set finalpsf $outName.psf
    set finalpdb $outName.pdb

    set pi [expr 4.0*atan(1.0)]
    set gamma [expr $angle*$pi/180.0]
    set m [expr tan($gamma)]

    set a [expr {5.0*$diamX}]
    set b [expr {5.0*$diamY}]
    #set selText "sqrt((x/($a + abs(z)*$m))^2 + (y/($b + abs(z)*$m))^2) < 1"

    set z0 [expr {5.0*$parabLen}]
    set a1 [expr {$a - 0.5*$m*$z0}]
    set b1 [expr {$b - 0.5*$m*$z0}]
    set g [expr {0.5*$m/$z0}]

    set selTextBig "sqrt((x/($a1 + abs(z)*$m))^2 + (y/($b1 + abs(z)*$m))^2) < 1"
    set selTextSmall "sqrt((x/($a + $g*z^2))^2 + (y/($b + $g*z^2))^2) < 1"
    set selText "(abs(z) > z0 and $selTextBig) or (abs(z) < $z0 and $selTextSmall)"
    #set selText $selTextBig

    # Obtain the {segid resid name} for the selection.
    mol load psf $psf pdb $pdb
    set sel [atomselect top $selText]
    foreach null {0} {
	set atomList [lsort -unique [$sel get {segname resid name}]]
	set posList [$sel get {x y z}]
    }
    set nAtoms [$sel num]
    set nResidues [llength $atomList]
    $sel delete

    package require psfgen 1.3
    resetpsf

    readpsf $psf
    coordpdb $pdb

    # Delete the selection.
    foreach atom $atomList {
	delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
    }

    writepsf tmp.psf
    writepdb tmp.pdb
    puts ""
    puts "$nAtoms atoms were deleted."
    puts "$nResidues residues were deleted."
    mol delete top


    # Find the largest radius cut.
    foreach {x y z} [lindex $posList 0] { break }
    set sMax [expr {sqrt($x*$x + $y*$y)}]
    foreach pos $posList {
	foreach {x y z} $pos { break }
	set s [expr {sqrt($x*$x + $y*$y)}]
	if {$s > $sMax} { set sMax $s }
    }
    puts "maximum radius: $sMax"

    set sMax2 [expr $sMax*$sMax]
    mol load psf tmp.psf pdb tmp.pdb
    set sel [atomselect top "($bufferText) and x^2+y^2<$sMax2"]
    foreach null {0} {
	set atomList [lsort -unique [$sel get {segname resid name}]]
    }
    puts "possible violators: [$sel num]"

    set violators {}
    foreach atom $atomList {
	if {rand() < $frac} {
	    lappend violators $atom
	}
    }
    set nAtoms [llength $violators]

    resetpsf
    readpsf tmp.psf
    coordpdb tmp.pdb

    # Delete the selection.
    foreach atom $violators {
	delatom [lindex $atom 0] [lindex $atom 1] [lindex $atom 2]
    }

    writepsf $finalpsf
    writepdb $finalpdb
    puts ""
    puts "$nAtoms atoms were deleted."
    mol delete top
}

exit
