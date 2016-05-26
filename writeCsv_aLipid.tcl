proc lavg L {expr ([join $L +])/[llength $L].}

proc RMSD {psfPrefix pdbPrefix dcd xstPrefix outPrefix} {

    variable selText
    variable stride
    variable timestep
    

    package require pbctools

    # Input:
    set psf $psfPrefix.psf
    set pdb $pdbPrefix.pdb

    # load trajectory and wait until loading completes 
    mol load psf $psf pdb $pdb

    set sel [atomselect top "$selText"]
    set all [atomselect top all]
    set lip [atomselect top "resname DFPC"]

    animate delete beg 0 end 1 skip 0 0

    mol addfile $dcd first 0 last -1 step $stride waitfor all
    set n [ molinfo top get numframes ]
    puts "  Loaded $n frames."
    
    # define first and last frames and stride
    set firstframe 0
    set lastframe $n
   
    # Open the output files.
    set out {}
    foreach dim {x y z} {
        lappend out [open "${outPrefix}.writeCsv_aLipid_${dim}.csv" w]
    }


    # Add index header
    foreach dim {x y z} o $out {

        puts -nonewline $o  " "

        set ind [$sel get index]

        foreach i $ind {
            puts -nonewline $o [format ",%d" ${i}]
        }
        puts $o ""

    }


    # Get the time in nanoseconds.
    set t ""
    set inStream [open ${xstPrefix}.xst r]
    gets $inStream;  gets $inStream;
    while {[gets $inStream line] > 0} {
        lappend t [expr [lindex $line 0] * $timestep * 1.0e-6]
    }


    # proceed frame by frame
    for {set f $firstframe} {$f <= $lastframe} {incr f} {
        
        $sel frame $f
        $all frame $f
        $lip frame $f

        # align dynamic selection to origin
        set alignVector [vecinvert [measure center $lip]]
        $all moveby $alignVector

        foreach dim {x y z} o $out {

            puts -nonewline $o  "[lindex $t $f]"

            set pos [$sel get $dim]

            foreach p $pos {
                puts -nonewline $o [format ",%.2f" ${p}]
            }
            puts $o ""

        }

        puts "frame $f" 
        
    }

    mol delete top

    foreach o $out {
        close $o
    }

}


if {$argc < 8} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcd xstPrefix selText stride outPrefix"
    exit
}

set psfPrefix    [lindex $argv 0]
set pdbPrefix    [lindex $argv 1]
set dcd          [lindex $argv 2]
set xstPrefix    [lindex $argv 3]

#map the '_' in selection text to ' '
set selText      [lindex $argv 4]
set selText      [string map {_ \ } $selText]

set stride       [lindex $argv 5]
set timestep     [lindex $argv 6]
set outPrefix    [lindex $argv 7]

RMSD $psfPrefix $pdbPrefix $dcd $xstPrefix $outPrefix 

exit



