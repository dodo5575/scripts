# Check that all of the files implied by trajectory index exist.
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc != 1} {
    puts "Usage: tclsh trajectoryExistence.tcl indexFile"
    exit
}


# Input:
set inFile [lindex $argv 0]

proc checkExistence {path} {
    if {![file exists $path]} {
	puts "WARNING! The file $path does not exist."
    }
}

set in [open $inFile r]
foreach line [split [read $in] \n] {
    # Ignore comments.
    if {[string match "\#*" $line]} {continue}
    
    set tok [concat $line]
    if {[llength $tok] < 9} {continue}

    foreach {index class name structDir structName xscName dcdDir dcdName dcdSet dcdFreq offset} $tok {break}

    set psf $structDir/$structName.psf
    checkExistence $psf
    set pdb $structDir/$structName.pdb
    checkExistence $pdb
    set xsc $structDir/$xscName.xsc
    checkExistence $xsc
    
    foreach d $dcdSet {
	set dcd $dcdDir/${dcdName}${d}.dcd
	checkExistence $dcd
    }
    if {[llength $dcdSet] == 0} {
	set dcd $dcdDir/${dcdName}.dcd
	checkExistence $dcd
    }
}

close $in



