#Program:
#(1)This program aligns the center of the molecule to the center of a specified box.
#(2)If the length of the box is set to be zero in all direction, it can align the system to the origin.
#Usage: 
#vmd -dispdev text -e $argv0 -args inpdbPrefix boxX boxY boxZ outpdbPrefix 
#Written by Chen-Yu Li    dodo5575@gmail.com
#2013/5/30


#puts "alignCenter {inpdbPrefix boxX boxY boxZ outpdbPrefix}"

proc vecSub {a b} {
    set c {}
    foreach ai $a bi $b {
        lappend c [expr $ai-$bi]
    }
    return $c
}

proc alignCenter {inpdbPrefix boxX boxY boxZ outpdbPrefix} {

	set boxVec [list $boxX $boxY $boxZ]

	set in [mol new ${inpdbPrefix}.pdb]

	set all [atomselect top all]	

	set centerOfMolecule [measure center $all]
	
	puts "The current center of the system is ${centerOfMolecule}."

	set centerOfBox {}

	foreach i $boxVec {
		lappend centerOfBox [expr $i / 2]

	}

	puts "The center of the specified box is ${centerOfBox}."	

	$all moveby [vecSub $centerOfBox $centerOfMolecule]
	
	if {$boxX!=0 && $boxY!=0 && $boxZ!=0} {
		molinfo top set a $boxX
		molinfo top set b $boxY
		molinfo top set c $boxZ
	}
	
	$all writepdb ${outpdbPrefix}.pdb

	mol delete $in
}

if {$argc < 5} {
	puts "vmd -dispdev text -e ~/scripts/alignCenter.tcl -args inpdbPrefix boxX boxY boxZ outpdbPrefix"
	exit
}

set inpdbPrefix [lindex $argv 0]
set boxX [lindex $argv 1]
set boxY [lindex $argv 2]
set boxZ [lindex $argv 3]
set outpdbPrefix [lindex $argv 4]

alignCenter $inpdbPrefix $boxX $boxY $boxZ $outpdbPrefix

puts "Alignment is complete."

exit
