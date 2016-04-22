# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname ADNA"
# Input:
set psf longDNA_AMBER.psf
set pdb longDNA_AMBER.pdb


mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

set resList [lsort -unique -integer [$sel get resid]]
$sel delete

set seqA ""
set seqB ""

foreach res $resList {
    set s [atomselect top "($selText) and resid $res"]
    set resName [lindex [$s get resname] 0]
    $s delete

    if {[string match "A*" $resName]} {
	set seqA ${seqA}A
	set seqB ${seqB}T
    } elseif {[string match "T*" $resName]} {
	set seqA ${seqA}T
	set seqB ${seqB}A
    } elseif {[string match "G*" $resName]} {
	set seqA ${seqA}G
	set seqB ${seqB}C
    } elseif {[string match "C*" $resName]} {
	set seqA ${seqA}C
	set seqB ${seqB}G
    } elseif {[string match "U*" $resName]} {
	set seqA ${seqA}U
	set seqB ${seqB}A
    } else {
	puts "Warning: Unknown residue ${resName}:${res}."
    }
}

puts "\nSequence: [llength $resList] nucleotides"
puts "sequence A:"
puts $seqA
puts "sequence B:"
puts [string reverse $seqB]

mol delete top
exit

