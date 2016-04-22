# Author: Jeff Comer <jcomer2@illinois.edu>

source $env(HOME)/scripts/useful.tcl
vmdargs getFit.tcl selText fromPdb toPdb outPdb outTransFile

set fromMol [mol new $fromPdb]
set fromSel [atomselect $fromMol $selText]

set toMol [mol new $toPdb]
set toSel [atomselect $toMol $selText]

set m [measure fit $fromSel $toSel weight mass]
$fromSel move $m
$fromSel writepdb trans.pdb
set m3 [matConvert4To3 $m]

puts ""
puts [lindex $m3 0 0]
puts [lindex $m3 0 1]
puts [lindex $m3 0 2]
puts [lindex $m3 1]
puts ""

set out [open $outTransFile w]
puts $out [lindex $m3 0 0]
puts $out [lindex $m3 0 1]
puts $out [lindex $m3 0 2]
puts $out [lindex $m3 1]
close $out

exit
