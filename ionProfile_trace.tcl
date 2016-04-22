proc getMinZ {xstPrefix} {

    set z 99999

    set ch [open $xstPrefix.xst] ;# open $file for reading
    while {[gets $ch line] > 0} { ## read each line of $ch
        ## $line = one of the following
        ## #$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w
        ## 0 40 0 0 0 40 0 0 0 40 0 0 0 0 0 0 0 0 0

        if {[string index $line 0] == "\#"} { continue } ;# skip comments

	if {[lindex $line 9] < $z} {
            set z [lindex $line 9]
	}
	
    }
    close $ch

    return $z
}

proc getAverageDimension {xstPrefix} {

    set xList {}
    set yList {}
    set zList {}

    set ch [open $xstPrefix.xst] ;# open $file for reading
    while {[gets $ch line] > 0} { ## read each line of $ch
        ## $line = one of the following
        ## #$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w
        ## 0 40 0 0 0 40 0 0 0 40 0 0 0 0 0 0 0 0 0

        if {[string index $line 0] == "\#"} { continue } ;# skip comments

        set x [lindex $line 1]
        set y [lindex $line 5]
        set z [lindex $line 9]

        if { [catch {set v [expr {$x*$y*$z}]} error] } {
            puts "x = $x"
            puts "y = $y"
            puts "z = $z"
            error $error
        }

        lappend xList $x
        lappend yList $y
        lappend zList $z
	
    }
    close $ch

    set length [llength $xList]

    set xTotal 0
    set yTotal 0
    set zTotal 0

    foreach x $xList y $yList z $zList {

	set xTotal [expr $xTotal + $x]
	set yTotal [expr $yTotal + $y]
	set zTotal [expr $zTotal + $z]

    }

    set xAve [expr $xTotal/$length]
    set yAve [expr $yTotal/$length]
    set zAve [expr $zTotal/$length]

    set systemSizeList [list $xAve $yAve $zAve]

    return $systemSizeList

}

#Compute the number of specific ion in a given section at a given frame.
proc computeIon {selIon sectionMin sectionMax} {

    set sel [atomselect top "$selIon and z >= $sectionMin and z < $sectionMax"]

    set num [$sel num]

    return $num
}


#Main function
proc compute {name structPrefix outDir dcdFreq timestep startFrame xstPrefix dcdList} {

    set displayPeriod 20
    set stride 1
    set selText [list "name POT" "name CLA" "name MG" "name P"]
    set nameList [list "K" "Cl" "Mg" "P"]

    set systemSize [getAverageDimension $xstPrefix]

    set xAve [lindex $systemSize 0]
    set yAve [lindex $systemSize 1]
    set zAve [lindex $systemSize 2]
    puts "$xAve $yAve $zAve"

    set MinZ [getMinZ $xstPrefix]

    set Zmax [expr $MinZ/2.0]
    set Zmin [expr -1*($MinZ/2.0)]

    #divide into section
    if {[expr floor($Zmax)] == $Zmax } {
	set floorMax [expr $Zmax - 1]
    } else {
	set floorMax [expr floor($Zmax)]
    }
    set floorMin [expr floor($Zmin+1)]
    set nSection [expr $floorMax - $floorMin + 2]

    # Input
    set psf $structPrefix.psf
    set pdb $structPrefix.pdb

    # Get the time change between frames in nanoseconds.
    set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

    # Open the output files.
    set out {}
    foreach n $nameList {
        lappend out [open "${outDir}/${name}_${n}_profile_trace.dat" w]
    }

    # Load the system.
    mol load psf $psf pdb $pdb

    set nuc [atomselect top "nucleic"]
    set all [atomselect top all]

    # Loop over the dcd files.
    set nFrames0 0
    foreach dcdFile $dcdList {
        # Load the trajectory.
        animate delete all
        mol addfile $dcdFile type dcd step $stride waitfor all
        set nFrames [molinfo top get numframes]
        puts [format "Reading %i frames." $nFrames]

        $nuc frame 0
        $all frame 0

        # Move forward, computing
        # current at each step.

	foreach s $selText o $out {
            for {set f $startFrame} {$f < $nFrames} {incr f} {

                if {$f == $startFrame} {
                    for {set secN 0} {$secN <= [expr $nSection - 1]} {incr secN} {

                        if {$secN == 0} {

                            set sectionMin $Zmin
                            set sectionMax $floorMin
	                    set interval [expr ($Zmin + $floorMin) / 2.0]

	                    set volume [expr ($sectionMax - $sectionMin) * $xAve * $yAve]

	                    puts -nonewline "$interval\t"
	                    puts -nonewline $o "$interval\t"
                            

                        } elseif {$secN == [expr $nSection - 1]} {
                     
                            set sectionMin $floorMax
                            set sectionMax $Zmax
	                    set interval [expr ($Zmax + $floorMax) / 2.0]

	                    set volume [expr ($sectionMax - $sectionMin) * $xAve * $yAve]

	                    puts -nonewline "$interval\n"
	                    puts -nonewline $o "$interval\n"

                        } else {

                            set sectionMax [expr $floorMin + $secN]
                            set sectionMin [expr $sectionMax - 1]
                            set interval [expr $sectionMax - 0.5]
	                    
	                    set volume [expr ($sectionMax - $sectionMin) * $xAve * $yAve]

	                    puts -nonewline "$interval\t"
	                    puts -nonewline $o "$interval\t"

                        }
	            }
                }

                for {set secN 0} {$secN <= [expr $nSection - 1]} {incr secN} {

                    if {$secN == 0} {

                	set sectionMin $Zmin
                	set sectionMax $floorMin
	                set interval [expr ($Zmin + $floorMin) / 2.0]

	                set volume [expr ($sectionMax - $sectionMin) * $xAve * $yAve]

                        animate goto $f

                        $nuc frame $f
                        $all frame $f

                        # Compute the number of ion for each selection.           
                        # Write the number of carriers and the current for this selection.
                        set selIon [computeIon $s $sectionMin $sectionMax]

	                #1M = 0.000602 molecule/angstrom^3
	                set molar [expr $selIon / ($volume * 0.000602)] 

	                #puts -nonewline "$molar\t"
	                puts -nonewline $o "$molar\t"

                    } elseif {$secN == [expr $nSection - 1]} {
                 
                	set sectionMin $floorMax
                	set sectionMax $Zmax
	                set interval [expr ($Zmax + $floorMax) / 2.0]

	                set volume [expr ($sectionMax - $sectionMin) * $xAve * $yAve]

                        animate goto $f

                        $nuc frame $f
                        $all frame $f

                        # Compute the number of ion for each selection.           
                        # Write the number of carriers and the current for this selection.
                        set selIon [computeIon $s $sectionMin $sectionMax]

	                #1M = 0.000602 molecule/angstrom^3
	                set molar [expr $selIon / ($volume * 0.000602)] 

	                #puts -nonewline "$molar\n"
	                puts -nonewline $o "$molar\n"

                    } else {

                	set sectionMax [expr $floorMin + $secN]
                	set sectionMin [expr $sectionMax - 1]
                	set interval [expr $sectionMax - 0.5]
	                
	                set volume [expr ($sectionMax - $sectionMin) * $xAve * $yAve]

                        animate goto $f

                        $nuc frame $f
                        $all frame $f

                        # Compute the number of ion for each selection.           
                        # Write the number of carriers and the current for this selection.
                        set selIon [computeIon $s $sectionMin $sectionMax]

	                #1M = 0.000602 molecule/angstrom^3
	                set molar [expr $selIon / ($volume * 0.000602)] 

	                #puts -nonewline "$molar\t"
	                puts -nonewline $o "$molar\t"

                    }

	        }
            }
	}

        set nFrames0 [expr $nFrames+$nFrames0]
    }

    foreach o $out {
        close $o
    }

}

if {$argc < 5} {
    puts "vmd -dispdev text -e $argv0 -args name structPrefix outDir dcdFreq timestep startFrame xstPrefix dcdFile0 \[dcdFile1...\]"
    exit
}

set name [lindex $argv 0]
set structPrefix [lindex $argv 1]
set outDir [lindex $argv 2]
set dcdFreq [lindex $argv 3]
set timestep [lindex $argv 4]
set startFrame [lindex $argv 5]
set xstPrefix [lindex $argv 6]
set dcdList [lrange $argv 7 end]

compute $name $structPrefix $outDir $dcdFreq $timestep $startFrame $xstPrefix $dcdList
exit

