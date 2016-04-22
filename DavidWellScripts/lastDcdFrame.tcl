source $env(SCRIPTS)/Procs.tcl

if { $argc < 4 } {
    puts stderr "Usage: vmd -dispdev text -e lastDcdFrame.tcl -args [ul psf] [ul dcd] [ul format] \[[ul frame]\]"
    exit -1
}

set psf [lindex $argv 0]
set dcd	[lindex $argv 1]
set pdb [string match pdb [lindex $argv 2]]

if { $argc == 5 } {
    set frame [lindex $argv 3]
    set out [file rootname $dcd]_frame$frame
}

set nframes [dcdFrames $dcd]
if { ! [info exists frame] } {
    set frame [expr $nframes-1]
}

set out [file rootname $dcd]_frame$frame
if { $pdb } {
    set out $out.pdb
} else  {
    set out $out.coor
}

mol load psf $psf
mol addfile $dcd first $frame last $frame

set all [atomselect top all]
if { ! [file exists $out] } {
    if { $pdb } {
	$all writepdb $out
    } else {
	$all writenamdbin $out
    }
} else {
    puts "$out already exists!"
}

exit
