source $env(SCRIPTS)/Procs.tcl

### ARGUMENTS ###
#
if { $argc < 4 || $argc > 6 } {
    puts stderr "Usage: vmd -dispdev text -e rmsd_per_res.tcl -args [ul psf] [ul dcd] [ul outbase] \[ \[ [ul first_n_frames] | [ul -last_n_frames] \] | [ul first_frame] [ul last_frame] \]"
    exit -1
}

set psf [shift argv]
set dcd [shift argv]
set outbase [shift argv]

set nframes_dcd [dcdFrames $dcd]

if { $argc == 4 } {
    puts "READING ALL FRAMES"
    set first 0
    set last [expr $nframes_dcd - 1]
} elseif { $argc == 5 } {
    set frame [shift argv]
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
    set first [shift argv]
    set last [shift argv]
}


### MAIN ###
#
set nframes [expr $last - $first + 1]
puts "READING FRAMES $first TO $last"

set mol [mol load psf $psf]
mol addfile $dcd first $first last $last waitfor all

# get list of segnames
set prot [atomselect $mol protein]
set segnames [lsort -unique [$prot get segname]]

foreach segname $segnames {
    puts "\nSEGNAME $segname"
    
    set out ${outbase}_segname${segname}.dat
    if { [file exists $out] } {
	puts "FILE $out ALREADY EXISTS, SKIPPING ..."
	continue
    } else {
	puts "OPENING FILE $out FOR WRITING"
    }
    set outch [open $out w]
    
    set seltext "protein name CA and segname $segname"
    set ref [atomselect $mol $seltext frame 0]
    set sel [atomselect $mol $seltext]
    
    # Align
    puts "ALIGNING"
    for { set i 0 } { $i < $nframes } { incr i } {
	progressbar $i $nframes 50
	$sel frame $i
	set mat [measure fit $sel $ref]
	$sel move $mat
    }
    
    set resids [lsort -unique -integer [$sel get resid]]
    
    puts "COMPUTING RMSF"
    set i 0
    set imax [llength $resids]
    foreach resid $resids {
	progressbar $i $imax 50
	incr i
	
	set sel2 [atomselect $mol "$seltext and resid $resid"]
	set rmsf [measure rmsf $sel2]
	
	puts $outch "$resid $rmsf"
	
	$sel2 delete
    }
    
    puts "CLOSING FILE $out"
    close $outch
    $sel delete
    $ref delete
}

exit
