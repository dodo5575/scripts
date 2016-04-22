# Check that all of the files implied by trajectory index exist.
# Author: Jeff Comer <jcomer2@illinois.edu>

set stride 1
set timestep 1.0
set displayPeriod 20
# Input:
set indexFile pore2.0_update_index.txt
#set indexSet {0 1 2}
set analysisScript stretchTraj.tcl
# Defines: stride timestep dcdFreq psf pdb xsc dcd outSuffix moiety displayPeriod

set inIndex [open $indexFile r]
foreach line [split [read $inIndex] \n] {
    # Ignore comments.
    if {[string match "\#*" $line]} {continue}
    
    set tok [concat $line]
    if {[llength $tok] < 9} {continue}

    foreach {index class name structDir structName xscName dcdDir dcdName dcdSet dcdFreq offset} $tok {break}

    # Locate the structure files.
    set psf $structDir/$structName.psf
    set pdb $structDir/$structName.pdb
    set xsc $structDir/$xscName.xsc

    # Locate the dcd files.
    set dcd {}
    foreach d $dcdSet {
	lappend dcd $dcdDir/${dcdName}${d}.dcd
    }
    if {[llength $dcdSet] == 0} {
	lappend dcd $dcdDir/${dcdName}.dcd
    }
    
    # Set the name of the output file.
    set outSuffix $name[lindex $dcdSet 0]-[lindex $dcdSet end].dat
    
    # Determine the type.
    if {[string match "*coil*" $class] || [string match "*loop*" $class]} {
	set moiety hpDNA
    } else {
	set moiety dsDNA
    }

    #if {[lsearch $indexSet $index] > 0} {
    puts "\n*****Analyzing $name with $analysisScript...\n"
    source $analysisScript
    #}
}
close $inIndex

puts "Analysis complete."
exit



