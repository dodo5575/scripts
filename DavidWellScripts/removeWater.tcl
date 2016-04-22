# extractDcd_protDnaIons.tcl
# author: dbwells2@uiuc.edu

source $env(SCRIPTS)/Procs.tcl


### PARAMETERS ###
#
set seltext	"not water"
set buffer	1000


### ARGUMENTS ###
#
vmdargs removeWater.tcl psf dcd


### MAIN ###
#
regsub "\.dcd$" $dcd "_nowater.dcd" outdcd_tmp; regsub "^output/" $outdcd_tmp "output-scratch/" outdcd
regsub "\.psf$" $psf "_nowater.psf" outpsf

if { [file exists $outdcd] } {
    puts "\nFILE $outdcd EXISTS, SKIPPING ... \n"
    exit
}

mol load psf $psf
set nframes [dcdFrames $dcd]
for { set i 0 } { $i * $buffer < $nframes } { incr i } {
    animate delete all
    set first [expr $i * $buffer]
    if { ($i + 1) * $buffer <= $nframes } {
	set last [expr ($i + 1) * $buffer - 1]
    } else {
	set last [expr $nframes - 1]
    }
    puts "READING FRAMES $first TO $last (TOTAL $nframes)"
    mol addfile $dcd first $first last $last waitfor all
    set sel [atomselect top $seltext]
    animate write dcd ${outdcd}_$i sel $sel
    lappend tmp_dcds ${outdcd}_$i
}

eval [concat catdcd -o $outdcd $tmp_dcds]

foreach tmp_dcd $tmp_dcds {
    puts "DELETING TEMPORARY FILE $tmp_dcd"
    file delete $tmp_dcd
}

if { ! [file exists $outpsf] } {
    puts "WRITING PSF FILE $outpsf ... \n"
    $sel writepsf $outpsf
}

exit
