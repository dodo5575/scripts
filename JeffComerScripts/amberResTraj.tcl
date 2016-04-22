# Calculate the distance between two CoMs for a trajectory.
# Writes step time(ns) position (nm) cross-sectional radius (nm) to a file.
# to use: vmd -dispdev text -e stretchDcd.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

# Requires: stride timestep dcdFreq psf pdb xsc dcd outSuffix moiety displayPeriod

switch $moiety {
    hpDNA {
	set selText "segname HAIR and resid 103 to 111 118 to 127"
	#set selText "segname HAIR and resid 130 140 150 160 170"
    }
    dsDNA {
	set selText "(segname ADNA and resid 5 to 54) or (segname BDNA and resid 107 to 58)"
    }
    default {
	puts stderr "ERROR: Unrecognized moiety $moiety."
	exit
    }
}

set alphaName {"O3'" "P" "O5'" "C5'"}
set alphaRes {-1 0 0 0}
set gammaName {"O5'" "C5'" "C4'" "C3'"}
set gammaRes {0 0 0 0}
set pi [expr 4.0*atan(1.0)]
set dt [expr 1.0e-6*$timestep*$dcdFreq*$stride]

# Load the system.
mol load psf $psf pdb $pdb

# Get the residue ids.
set sel [atomselect top "$selText"]
set residue [lsort -unique [$sel get {segname resid}]]
$sel delete
set alpha {}
set gamma {}
set goodRes {}
foreach res $residue {
    set segName [lindex $res 0]
    set resId [lindex $res 1]
    
    set aList {}
    set gList {}
    foreach a $alphaName g $gammaName ar $alphaRes gr $gammaRes {
    	set sa [atomselect top "segname $segName and resid [expr $resId+$ar] and name $a"]
	if {[$sa num] == 1} {
	    lappend aList [lindex [$sa get index] 0]
	} else {
	    puts "Error: Atom ($segName $resId $a) does not exist."
	}
	$sa delete

	set sg [atomselect top "segname $segName and resid [expr $resId+$gr] and name $g"]
	if {[$sg num] == 1} {
	    lappend gList [lindex [$sg get index] 0]
	} else {
	    puts "Error: Atom ($segName $resId $g) does not exist."
	}
	$sg delete
    }

    lappend goodRes $resId
    lappend alpha $aList
    lappend gamma $gList
}

# Open the output files.
set outAlpha {}
set outGamma {}
foreach r $goodRes {
    lappend outAlpha [open amber/alpha_res${r}_${outSuffix} w]
    lappend outGamma [open amber/gamma_res${r}_${outSuffix} w]
}

# Loop over the dcd files.
set nFrames0 0
foreach dcdFile $dcd {
    # Load the trajectory.
    animate delete all
    mol addfile $dcdFile type dcd step $stride waitfor all
    set nFrames [molinfo top get numframes]
    puts [format "Reading %i frames." $nFrames]

    # Move forward computing the center-of-mass at every step.
    for {set f 0} {$f < $nFrames} {incr f} {
	molinfo top set frame $f

	# Get the time in nanoseconds for this frame.
	set t [expr ($nFrames0+$f)*$dt]

	# Determine the torsions.
	foreach a $alpha oa $outAlpha {
	    set alpha0 [measure dihed $a]
	    puts $oa "$t $alpha0"
	}
	foreach g $gamma og $outGamma {
	    set gamma0 [measure dihed $g]
	    puts $og "$t $gamma0"
	}
		
	# Write the time and distance.
	if {$f % $displayPeriod == 0} {
	    puts -nonewline [format "FRAME %i: " $f]
	    puts "$t $alpha0 $gamma0"
	}
    }
    set nFrames0 [expr $nFrames+$nFrames0]
}

foreach oa $outAlpha og $outGamma {
    close $oa
    close $og
}

mol delete top



