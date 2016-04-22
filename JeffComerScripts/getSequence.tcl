# Get the system size using measure minmax.
# to use: vmd -dispdev text -e getSize.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname PRTA PRTB"
# Input:
set psf BamHI_nonspec-0.1M.ver2a.psf
set pdb BamHI_nonspec-0.1M.ver2a.pdb


mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

set resList [lsort -unique -integer [$sel get resid]]
$sel delete

set seq ""
set nucleo {{ADE A} {THY T} {GUA G} {CYT C} {URA U}}
set aminoBasic {{LYS K} {ARG R} {HIS H} {HSE H} {HSD H} {HSP H}}
set aminoAcidic {{ASP D} {GLU E}}
set aminoPolar {{ASN N} {GLN Q} {SER S} {THR T} {TYR Y}}
set aminoNonpolar {{ALA A} {VAL V} {LEU L} {ILE I} {PRO P} {PHE F} {MET M} {TRP W} {GLY G} {CYS C}}
set nameList [concat $nucleo $aminoBasic $aminoAcidic $aminoPolar $aminoNonpolar]

foreach res $resList {
    set s [atomselect top "($selText) and resid $res"]
    set resName [lindex [$s get resname] 0]
    $s delete
    
    set unrec 1
    foreach name $nameList {
	if {[string match $resName [lindex $name 0]]} {
	    set seq "${seq}[lindex $name 1]"
	    set unrec 0
	    break
	}
    }
    
    if {$unrec} {
	set seq "${seq}?"
	puts "Warning: Unknown residue ${resName}:${res}."
    }
}

puts "\nSequence length: [llength $resList]"
puts "Sequence:"
puts $seq

mol delete top
exit



