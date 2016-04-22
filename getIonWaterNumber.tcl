# Calculate the number of ions, phosphorus and water
# Usage: vmd -dispdev text -e $argv0 -args name structPrefix
# Chen-Yu Li     cli56@illinois.edu
# 2014/1/27


if {$argc < 2} {
    puts "vmd -dispdev text -e $argv0 -args name structPrefix"
    exit
}

set name [lindex $argv 0]
set structPrefix [lindex $argv 1]

mol load psf ${structPrefix}.psf pdb ${structPrefix}.pdb

set out [open ${name}_ionNum.dat w]

set P [atomselect top "name P"]

set Cl [atomselect top "name CLA"]

set K [atomselect top "name POT"]

set Mg [atomselect top "name MG"]

set WT [atomselect top "name OH2"]

set PN [$P num]

set ClN [$Cl num]

set KN [$K num]

set MgN [$Mg num]

set WTN [$WT num]

puts $out "P\t$PN" 
puts $out "Cl\t$ClN" 
puts $out "K\t$KN" 
puts $out "Mg\t$MgN" 
puts $out "water\t$WTN"

close $out

exit
