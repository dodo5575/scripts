proc lavg L {expr ([join $L +])/[llength $L].}

proc RMSD {psfPrefix pdbPrefix dcd xstPrefix lipidResname dim R timestep outPut} {

    package require pbctools

    # Input:
    set psf $psfPrefix.psf
    set pdb $pdbPrefix.pdb

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
    mol load psf $psf pdb $pdb

    set lipidN [atomselect top "resname $lipidResname and name N"]
    $lipidN frame 0

    set lipidC216 [atomselect top "resname $lipidResname and name C216"]
    $lipidC216 frame 0
    set lipidC316 [atomselect top "resname $lipidResname and name C316"]
    $lipidC316 frame 0

    set C216_z_list [$lipidC216 get $dim]
    set C316_z_list [$lipidC316 get $dim]

    set residue_list [$lipidN get residue]

    set tail [vecscale [vecadd $C216_z_list $C316_z_list] 0.5]

    set vecRef [vecsub [$lipidN get $dim] $tail]

    animate delete beg 0 end 1 skip 0 0

    mol addfile $dcd first 0 last -1 waitfor all
    set n [ molinfo top get numframes ]
    puts "  Loaded $n frames."
    
    # define first and last frames and stride
    set firstframe 0
    set lastframe $n
    set stride 1
   
    # set output file
    set out [open $outPut w]

    # Get the time in nanoseconds.
    set t ""
    set inStream [open ${xstPrefix}.xst r]
    gets $inStream;  gets $inStream;
    while {[gets $inStream line] > 0} {
        lappend t [expr [lindex $line 0] * $timestep * 1.0e-6]
    }


    # proceed frame by frame
    for {set f $firstframe} {$f <= $lastframe} {set f [ expr $f + $stride ] } {
        
        $lipidN frame $f
        $lipidC216 frame $f
        $lipidC316 frame $f

        set C216_z_list [$lipidC216 get $dim]
        set C316_z_list [$lipidC316 get $dim]

        set tail [vecscale [vecadd $C216_z_list $C316_z_list] 0.5]

        set vecCurrent [vecsub [$lipidN get $dim]  $tail]

        set mul [vecmul $vecRef $vecCurrent]

        set resi_list ""
        set count 0
        foreach i $mul r $residue_list {
            if {$i < 0} {
                incr count
                lappend resi_list $r
            }
        }

        set sel [atomselect top "resname DFPC and same residue as within ${R} of nucleic"]
        set num [llength [lsort -unique [$sel get residue]]]
        $sel delete

        puts $out "[lindex $t $f]\t$count\t$num\t$resi_list"
        
        puts "frame $f" 
        
    }

    mol delete top

    close $out

}


if {$argc < 8} {
    puts "vmd -dispdev text -e $argv0 -args psfPrefix pdbPrefix dcd xstPrefix lipidResname dim R timestep outPut"
    exit
}

set psfPrefix    [lindex $argv 0]
set pdbPrefix    [lindex $argv 1]
set dcd          [lindex $argv 2]
set xstPrefix    [lindex $argv 3]
set lipidResname [lindex $argv 4]
set dim          [lindex $argv 5]
set R            [lindex $argv 6]
set timestep     [lindex $argv 7]
set outPut       [lindex $argv 8]

RMSD $psfPrefix $pdbPrefix $dcd $xstPrefix $lipidResname $dim $R $timestep $outPut

exit



