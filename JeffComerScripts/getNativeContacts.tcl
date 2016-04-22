# Extract native contacts.
# To run: vmd -dispdev text -e getNativeContacts.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

set selText "segname HAIR"
set dCut 9.0
set exclude {{SUG PHO 0} {SUG PHO 1} {SUG AB 0} {SUG TB 0} {SUG CB 0} {SUG GB 0}}
set dResMax 10000
set dResMin 0

#Input:
set psf cg_hairpin_neutral.psf
set coords cg_hairpin_native.pdb
#Output:
set outFile hairpin_native_contacts.txt

mol load psf $psf
mol addfile $coords
set sel [atomselect top $selText]
foreach zero {0} {
    set pos [$sel get {x y z}]
    set atom [$sel get {segname resid name}]
    set index [$sel get index]
}
set n [$sel num]
$sel delete

set out [open $outFile w]
# Find all native contacts.
set nm1 [expr $n - 1]
set nContacts 0
for {set i 0} {$i < $nm1} {incr i} {
    for {set j [expr $i+1]} {$j < $n} {incr j} {
	set iPos [lindex $pos $i]
	set jPos [lindex $pos $j]
	set d [veclength [vecsub $jPos $iPos]]

	if {$d < $dCut} {
	    set iId [lindex $atom $i]
	    set jId [lindex $atom $j]
	    set iName [lindex $iId 2]
	    set jName [lindex $jId 2]
	    set iRes [lindex $iId 1]
	    set jRes [lindex $jId 1]

	    # Check that this pair is not excluded.
	    set valid 1
	    foreach e $exclude {
		foreach {n0 n1 dRes} $e {break}

		if {$iRes + $dRes == $jRes} {
		    if {[string equal $iName $n0] && [string equal $jName $n1]} {
			set valid 0
			break
		    }
		}
		
		if {$jRes + $dRes == $iRes} {
		    if {[string equal $jName $n0] && [string equal $iName $n1]} {
			set valid 0
			break
		    }
		}
	    }
	    # Check the residue requirements.
	    if {abs($iRes - $jRes) > $dResMax} {set valid 0}
	    if {abs($iRes - $jRes) < $dResMin} {set valid 0}

	    if {$valid} {
		puts "$iId $jId $d"
		puts $out "$iId $jId $d"
		incr nContacts
	    }
	}
    }
}
close $out

puts "Wrote $nContacts native contacts to $outFile."
exit



