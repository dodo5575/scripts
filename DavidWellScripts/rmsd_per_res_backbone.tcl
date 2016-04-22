source $env(SCRIPTS)/Procs.tcl

### ARGUMENTS ###
#
if { $argc < 5 || $argc > 7 } {
    puts stderr "Usage: vmd -dispdev text -e rmsd_per_res.tcl -args [ul psf] [ul ref_coords] [ul dcd] [ul outbase] \[ \[ [ul first_n_frames] | [ul -last_n_frames] \] | [ul first_frame] [ul last_frame] \]"
    exit -1
}

set psf [lindex $argv 0]
set refpdb [lindex $argv 1]
set dcd [lindex $argv 2]
set outbase [lindex $argv 3]

set nframes_dcd [dcdFrames $dcd]

if { $argc == 5 } {
    puts "READING ALL FRAMES"
    set first 0
    set last [expr $nframes_dcd - 1]
} elseif { $argc == 6 } {
    set frame [lindex $argv 4]
    if { $frame > 0 } {
	# take first $frame frames
	puts "READING FIRST $frame FRAMES"
	set first 0
	set last [expr $frame - 1]
    } else {
	# take last $frame frames
	set frame [expr -$frame]
	puts "READING LAST $frame FRAMES"
	set first [expr $nframes_dcd - $frame]
	set last [expr $nframes_dcd - 1]
    }
} else {
    set first [lindex $argv 4]
    set last [lindex $argv 5]
}


### MAIN ###
#
set nframes [expr $last - $first + 1]
puts "READING FRAMES $first TO $last"

set mol [mol load psf $psf]
mol addfile $dcd first $first last $last waitfor all

set refmol [mol load psf $psf]
mol addfile $refpdb waitfor all

# get list of segnames
set prot [atomselect $mol protein]
set segnames [lsort -unique [$prot get segname]]

foreach segname $segnames {
    puts "\nSEGNAME $segname"
    
    set out ${outbase}_segname${segname}.dat
    if { [file exists $out] } {
	puts "FILE $out ALREADY EXISTS, SKIPPING ..."
	continue
    }
    set outch [open $out w]
    
    set seltext "protein backbone and segname $segname"
    set ref [atomselect $refmol $seltext]
    set sel [atomselect $mol $seltext]
    
    puts $outch "# psf: $psf"
    puts $outch "# refpdb: $refpdb"
    puts $outch "# dcd: $dcd"
    puts $outch "# first frame: $first"
    puts $outch "# last frame: $last"
    puts $outch "# seltext: $seltext"

    trajdev2 $ref $sel
    
    set resids [lsort -unique -integer [$sel get resid]]
    foreach resid $resids {
	set sel2 [atomselect $mol "$seltext and resid $resid"]
	set dev2s {}
	
	for { set frame 0 } { $frame < $nframes } { incr frame } {
	    $sel2 frame $frame
	    lappend dev2s [$sel2 get user]
	}
	
	set dev2s [eval concat $dev2s]
	set meandev2 [lmean $dev2s]
	set err [lstddev $dev2s]
	
	#puts "$resid $mean $err"
	puts $outch "$resid [expr sqrt($meandev2)] [expr 0.5 * $err / sqrt($meandev2)]"
	
	$sel2 delete
    }
    
    close $outch
    $sel delete
    $ref delete
}

exit
