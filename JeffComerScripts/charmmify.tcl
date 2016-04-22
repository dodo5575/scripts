# to use: vmd -dispdev text -e typeName.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set nameChange {{"resname THY and name C7" C5M} {"resname THY and name H71" H51} {"resname THY and name H72" H52} {"resname THY and name H73" H53}}
set selText "segname ADNA BDNA"
# Input:
set psf str_dna_int2.psf
set pdb str_dna_int2.pdb
set nameTypeFile name-type_charmm.txt
# Output:
set outPrefix str_dna_charmmed

# Get the name->type conversion.
set in [open $nameTypeFile r]
set conName {}
set conType {}
set conCharge {}
foreach line [split [read $in] \n] {
    if {[string length $line] < 2} {continue}
    if {[string match "\#*" $line]} {continue}
    set tok [concat $line]
    lappend conName [lrange $tok 0 1]
    lappend conType [lindex $tok 2]
    lappend conCharge [lindex $tok 3]
}
close $in

# Load the system.
mol load psf $psf pdb $pdb
set sel [atomselect top $selText]

# Change the atom names.
foreach change $nameChange {
    foreach {text name} $change {break}
    set s [atomselect top $text]
    $s set name $name
    puts "Changed the name of [$s num] atoms defined by ($text) to $name."
    $s delete
}

# Change the atom types using the atom names.
foreach zero {0} {set nameList [$sel get {resname name}]}
set newType {}
set newCharge {}
foreach n $nameList {
    set i [lsearch $conName $n]
    if {$i < 0} {
	puts "Warning! No corresponding type found for name $n."
	continue
    }
    lappend newType [lindex $conType $i]
    lappend newCharge [lindex $conCharge $i]
}
$sel set type $newType
$sel set charge $newCharge
$sel delete

set all [atomselect top all]
$all writepsf $outPrefix.psf
$all writepdb $outPrefix.pdb
$all delete

mol delete top
exit

