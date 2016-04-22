

proc RMSD {psfPrefix dcdPrefix xstPrefix lipidResname dim cutOff D timestep outPut} {

    package require pbctools

    # Input:
    set psf $psfPrefix.psf
    set dcd $dcdPrefix.dcd

    switch -- $dim {
    
        "x" {
            set dimension 0
            set a "y"
            set b "z"
        }
        "y" {
            set dimension 1
            set a "x"
            set b "z"
        }
        "z" {
            set dimension 2
            set a "x"
            set b "y"
        }
    }

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
        
        #set lipid [atomselect top "resname $lipidResname and same residue as within $D of nucleic" frame $i]

        # this is faster
        set lipid [atomselect top "resname $lipidResname and name P and ($a*$a + $b*$b > $D*$D)" frame $i] 

        set residue [$lipid get residue]

        #puts "$resID"

        set count 0

        foreach id $residue {

            set sel1 [atomselect top "resname $lipidResname and residue $id and name P" frame $i]
            set sel2 [atomselect top "resname $lipidResname and residue $id and name C216 C316" frame $i]

            set c1 [measure center $sel1]
            set c2 [measure center $sel2]

            set x1 [lindex $c1 0]
            set y1 [lindex $c1 1]
            set x2 [lindex $c2 0]
            set y2 [lindex $c2 1]

            if {[expr ($x2*$x2 + $y2*$y2) - ($x1*$x1 + $y1*$y1)] > [expr $cutOff*$cutOff]} {
                incr count
            } 
            $sel1 delete
            $sel2 delete

        }

        set percentage [expr $count / double([llength $residue])]

        puts $out "[lindex $t $i]\t$percentage"
        
        puts "frame $i" 
        #puts "frame $i $mid $up $down $a $b $D [$sel num]" 
        $lipid delete
        
    }

    mol delete top

    close $out

}


if {$argc < 8} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix dcdPrefix xstPrefix lipidResname dim cutOff D timestep outPut"
    exit
}

set psfPrefix    [lindex $argv 0]
set dcdPrefix    [lindex $argv 1]
set xstPrefix    [lindex $argv 2]
set lipidResname [lindex $argv 3]
set dim          [lindex $argv 4]
set cutOff       [lindex $argv 5]
set D            [lindex $argv 6]
set timestep     [lindex $argv 7]
set outPut       [lindex $argv 8]

RMSD $psfPrefix $dcdPrefix $xstPrefix $lipidResname $dim $cutOff $D $timestep $outPut

exit



