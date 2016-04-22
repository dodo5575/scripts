set x0 -4
set x1 2
set y0 -1
set y1 4
set dz 1.89355184
set startZ -11
set endZ 11
set instPrefix instance/nuc1_cyt
set instSuffix dat
# Input:
set prefix comb
set runScript ./meanAutocorrVelXY.tcl
# Output:
set outFile segmentAutocorr.sh
set corrPrefix correlation/autocorrxy
set indexFile segment_index.txt

mol load psf $prefix.psf pdb $prefix.pdb
set out [open $outFile w]
set outInd [open $indexFile w]
set index 0

foreach ion {pot chl} {
    for {set z $startZ} {$z <= $endZ} {set z [expr {$z + $dz}]} {
	set z0 [expr {$z - 0.5*$dz}]
	set z1 [expr {$z + 0.5*$dz}]
	set selText "segname NULL and x > $x0 and x < $x1 and y > $y0 and y < $y1 and z > $z0 and z < $z1"
	set sel [atomselect top "$selText"]

	# Get the node list.
	foreach quiet {0} {set nodeList [$sel get resid]}

	# Write the index data.
	puts "\nlevel $index $z"
	puts $nodeList
	puts $outInd "\nlevel $index $z"
	puts $outInd $nodeList

	# Write the script.
	set fileList {}
	foreach node $nodeList {
	    lappend fileList ${instPrefix}_${ion}.$node.$instSuffix
	}
	puts $out "${runScript} ${corrPrefix}_${ion}.$index $fileList"

	$sel delete
	incr index
    }
}

close $out
close $outInd
mol delete top
exit
