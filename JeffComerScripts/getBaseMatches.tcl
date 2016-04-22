# Extract native contacts.
# To run: vmd -dispdev text -e getBaseMatches.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname HAIR"
set match {{AB TB} {CB GB}}

#Input:
set psf cg_hairpin_neutral.psf
set coords cg_hairpin_neutral.pdb
#Output:
set outFile hairpin_base_matches.txt

mol load psf $psf
mol addfile $coords

set nMatches 0
set out [open $outFile w]
# Find all possible matching bases.
foreach m $match {
    set selA [atomselect top "($selText) and name [lindex $m 0]"]
    set selB [atomselect top "($selText) and name [lindex $m 1]"]
    set atomA [$selA get {segname resid name}]
    set atomB [$selB get {segname resid name}]
    set nMatches [expr $nMatches + [$selA num]*[$selB num]]
    $selA delete
    $selB delete
    
    foreach a $atomA {
	foreach b $atomB {
	    puts $out "$a $b"
	}
    }
} 
close $out

puts "Wrote $nMatches base matches to $outFile."
exit



