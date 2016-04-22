#Attention: use '_' instead of ' ' in selection text

proc RMSD {psfPrefix pdbPrefix dcdPrefix xstPrefix selText timestep outPut} {

    package require pbctools

    # Input:
    set psf $psfPrefix.psf
    set pdb $pdbPrefix.pdb
    set dcd $dcdPrefix.dcd

    # set reference selection
    mol load psf $psf pdb $pdb
    set selref [atomselect top $selText]
    puts "number of atoms in the selection: [$selref num]"

    # load trajectory and wait until loading completes 
    mol load psf $psf dcd $dcd
    set n [ molinfo top get numframes ]
    puts "  Loaded $n frames."
    
    # define first and last frames and stride
    set firstframe 0
    set lastframe $n
    set stride 1
   
    # set output file
    set out [open $outPut a]

    # Get the time in nanoseconds.
    set t ""
    set inStream [open ${xstPrefix}.xst r]
    gets $inStream;  gets $inStream;
    while {[gets $inStream line] > 0} {
        lappend t [expr [lindex $line 0] * $timestep * 1.0e-6]
    }
 
    # proceed frame by frame
    for {set i $firstframe} {$i <= $lastframe} {set i [ expr $i + $stride ] } {
        
        set seldyn [atomselect top $selText frame $i]

        puts $out "[lindex $t $i]\t[measure center $seldyn ]"
        
        $seldyn delete
        
        puts "frame $i" 
    }

    mol delete top

    close $out

}


if {$argc < 7} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcdPrefix xstPrefix selText timestep outPut"
    exit
}

set psfPrefix [lindex $argv 0]
set pdbPrefix [lindex $argv 1]
set dcdPrefix [lindex $argv 2]
set xstPrefix [lindex $argv 3]
#map the '_' in selection text to ' '
set selText [lindex $argv 4]
set selText   [string map {_ \ } $selText]
set timestep [lindex $argv 5]
set outPut [lindex $argv 6]

RMSD $psfPrefix $pdbPrefix $dcdPrefix $xstPrefix $selText $timestep $outPut

exit



