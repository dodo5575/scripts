# Author: jcomer2@uiuc.edu

set selText "resname SIN"
set radiusList {"type \"SI.*\"" 2.135 "type \"N.*\"" 1.9975}
set buffer 3.0
# Input:
set psf open2.0_0.1M.psf
set pdb open_trap2.0_1V2.pdb
set xsc open_trap_1V0.restart.xsc
# Output:
set outFile open_density.dx

mol load psf $psf pdb $pdb

# Set the radii.
foreach {tex rad} $radiusList {
    set sel [atomselect top "($selText) and ($tex)"]
    $sel set radius $rad
    $sel delete
}


# Read the system size from the xsc file.
set in [open $xsc r]
foreach line [split [read $in] "\n"] {
    if {![string match "#*" $line]} {
	set param [split $line]
	puts $param
	set lx [lrange $param 1 3]
        set ly [lrange $param 4 6]
	set lz [lrange $param 7 9]
	break
    }
}
puts "NOTE: The system size is ($lx) ($ly) ($lz).\n"
close $in

# Create the volume map.
set origin [vecadd [vecscale -0.5 $lx] [vecscale -0.5 $ly] [vecscale -0.5 $lz]]
set destin [vecadd [vecscale 0.5 $lx] [vecscale 0.5 $ly] [vecscale 0.5 $lz]]
set sel [atomselect top $selText]
set minmax [measure minmax $sel]
set origin [vecsub [lindex $minmax 0] [list $buffer $buffer $buffer]]
set destin [vecadd [lindex $minmax 1] [list $buffer $buffer $buffer]]
volmap density $sel -minmax [list $origin $destin] -checkpoint 0 -o $outFile -allframes -combine avg
exit
