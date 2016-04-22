#!/usr/bin/tclsh
# Split a file of mini-trajectories into many files differing by position in space. 
# Author: Jeff Comer <jcomer2@illinois.edu>

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

if {$argc < 6} {
    puts "Usage: ./assignInstancesToGridBinary.tcl skipSteps instanceSteps dt gridFile outDir outName inFileList"
    exit
}

set skipSteps [lindex $argv 0]
set instanceSteps [lindex $argv 1]
set dt [lindex $argv 2]
set gridFile [lindex $argv 3]
set outDir [lindex $argv 4]
set outName [lindex $argv 5]
set inFileList [lrange $argv 6 end]

if {![string is integer $skipSteps]} {
    puts "ERROR: The skipSteps must be an integer."
    exit
}
if {![string is integer $instanceSteps]} {
    puts "ERROR: The instanceSteps must be an integer."
    exit
}
if {![string is double $dt]} {
    puts "ERROR: The dt must be a double."
    exit
}

proc writeList {fileName l} {
    set out [open $fileName w]
    foreach i $l {
	puts $out $i
    }
    close $out
    return
}

proc writeCorrFunction {fileName l dt} {
    set out [open $fileName w]
    set step 0
    foreach i $l {
	puts $out "[expr {$dt*$step}] $i"
	incr step
    }
    close $out
    return
}

proc meanPos {posList} {
    if {[llength $posList] == 0} { return [vecZero] }

    set sum [vecZero]
    foreach pos $posList {
	set sum [vecAdd $sum $pos]
    }
    set sum [vecScale [expr 1.0/[llength $posList]] $sum]
    return $sum
}

proc wrapToSelf {gridVar posList} {
    upvar $gridVar grid

    if {[llength $posList] == 0} { return $posList }
    set r0 [lindex $posList 0]

    set ret {}
    foreach pos $posList {
	set dr [vecSub $pos $r0]
	set r1 [vecAdd $r0 [wrapDiff grid $dr]]
	lappend ret $r1
    }
    return $ret
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr {$ind-1}]]
}

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}


# Here we assume a timestep of 1 fs.
# This gives the interval between the positions in the binary output.
proc makeVelocityForce {gridVar pos0 pos1 force t} {
    upvar $gridVar grid
    
    set d [wrapDiff grid [vecSub $pos1 $pos0]]
    set line "$t $d $force"
    return $line
}

proc distributeCorr {gridVar inFileList outDir outName skipSteps trajSteps dt} {
    upvar $gridVar grid
    set dtInv [expr {1.0/$dt}]

    foreach inFile $inFileList {
	# Open the files.
	set in [open $inFile r]
	fconfigure $in -translation binary
	set chunk 36
	set inBinData [read $in]	
	puts "Read $inFile."

	set count 0
	set lineList {}
	set posList {}
	set readCount 0
	while {[string length $inBinData] >= $chunk} {
	    incr readCount
	    if {$readCount <= $skipSteps} { continue }

	    set fields [binary scan $inBinData f3f3f3 pos0 pos1 force]
	    set line [makeVelocityForce grid $pos0 $pos1 $force [expr {$dt*$readCount}]]
	    set inBinData [string range $inBinData $chunk end]

	    if {[llength $lineList] >= $trajSteps} {
		# Find the node of this instance.
		# Watch for the periodic boundaries.
		set cen [meanPos [wrapToSelf grid $posList]]
		set node [nearestIndex grid $cen]

		# Copy the instance into the file of its grid node.
		lappend lineList END
		set fileName "$outDir/$outName.$node.dat"
		set out [open $fileName a]
		foreach l $lineList {
		    puts $out $l
		}
		close $out
		
		set lineList {}
		set posList {}
		incr count
		continue
	    }
	    
	    # Build up the list.
	    lappend lineList $line
	    lappend posList $pos0
	}
	close $in
    }
    return $count
}

# Get the grid geometry.
readDx grid $gridFile
puts "Read the grid with dimensions $grid(nx) $grid(ny) $grid(nz)."

# Delete the files in the output directory that would collide.
set collideList [glob -nocomplain $outDir/$outName.*.dat]
if {[llength $collideList] > 0} {
    foreach collide $collideList {
	file delete $collide
    }
}

# Sequester each instance (mini-trajectory) in a file for the corresponding grid node.
set n [distributeCorr grid $inFileList $outDir $outName $skipSteps $instanceSteps $dt]
puts "Distributed $n instances."

exit
