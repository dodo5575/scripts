#!/usr/bin/tclsh
# Split a file of mini-trajectories into many files differing by position in space. 
# Author: Jeff Comer <jcomer2@illinois.edu>

source $env(HOME)/scripts/vector.tcl
source $env(HOME)/scripts/gridForce.tcl

if {$argc < 8} {
    puts "Usage: ./assignInstancesToCyl.tcl isBinary skipSteps instanceSteps dt cylinderFile systemGrid outDir outName logFile0 logFile1..."
    puts "isBinary: 0 for NAMD log file input, 1 for binary file input"
    puts "skipSteps is the number of records to skip before recording instances for each file."
    puts "instanceSteps are the number of records in the mini-trajectories that will be output."
    puts "dt is the time between records in the log files."
    puts "cylinderFile defines the cylindrical grid: "
    puts "s0 z0"
    puts "ds dz"
    puts "ns nz"
    puts "systemGrid is the dx file that gives the system size--used for wrapping."
    puts "outDir is the directory in which the output files go."
    puts "outName is the prefix of the output files."
    puts "logFile* are NAMD logFiles containing SMD lines."
    puts "If isBinary != 0 they are binary files containing velocity and force data created by `tclForcesData.tcl'."
    exit
}

set isBinary [lindex $argv 0]
set skipSteps [lindex $argv 1]
set instanceSteps [lindex $argv 2]
set dt [lindex $argv 3]
set cylFile [lindex $argv 4]
set gridFile [lindex $argv 5]
set outDir [lindex $argv 6]
set outName [lindex $argv 7]
set inFileList [lrange $argv 8 end]

if {![string is integer $isBinary]} {
    puts "ERROR: isBinary must be an integer."
    exit
}
if {![string is integer $skipSteps]} {
    puts "ERROR: skipSteps must be an integer."
    exit
}
if {![string is integer $instanceSteps]} {
    puts "ERROR: instanceSteps must be an integer."
    exit
}
if {![string is double $dt]} {
    puts "ERROR: dt must be a double."
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


# This gives the velocity between the positions.
proc makeVelocityForce {gridVar pos0 pos1 force t dt} {
    upvar $gridVar grid
    
    set d [wrapDiff grid [vecScale [expr {1.0/$dt}] [vecSub $pos1 $pos0]]]
    set line "$t $d $force"
    return $line
}

proc cylGeoRead {cylVar fileName} {
    upvar $cylVar cyl

    set cyl(s0) 0.0
    set cyl(z0) 0.0
    set cyl(ds) 1.0
    set cyl(dz) 1.0
    set cyl(ns) 1
    set cyl(nz) 1

    set list0 {s0 ds ns}
    set list1 {z0 dz nz}

    set count 0
    set in [open $fileName r]
    while {[gets $in line] >= 0} {
	if {[string match "\#*" $line]} { continue }

	set tok [concat $line]

	if {$count >= [llength $list0]} { continue }
	set cyl([lindex $list0 $count]) [lindex $tok 0]
	set cyl([lindex $list1 $count]) [lindex $tok 1]
	incr count
    }
    close $in
}

proc cylGeoNearestIndex {cylVar r} {
    upvar $cylVar cyl
    
    foreach {x y z} $r { break }
    set s [expr {sqrt($x*$x+$y*$y)}]
    set z [expr {abs($z)}]
    
    set is [expr {int(floor(($s-$cyl(s0))/$cyl(ds)))}]
    set iz [expr {int(floor(($z-$cyl(z0))/$cyl(dz)))}]

    if {$is < 0 || $is >= $cyl(ns)} { return -1 }
    if {$iz < 0 || $iz >= $cyl(nz)} { return -1 }
    return [expr {$iz + $cyl(nz)*$is}]
}

proc distributeCorrCyl {gridVar cylVar inFileList outDir outName skipSteps trajSteps dt} {
    upvar $gridVar grid
    upvar $cylVar cyl

    foreach inFile $inFileList {
	# Open the files.
	set in [open $inFile r]
	set inData [open $inFile r]
	puts "Reading `$inFile...'"
	
	set count 0
	set lineList {}
	set posList {}
	set readCount 0
	set pos0 [vecZero]
	while {[gets $in line] >= 0} {
	    if {![string match "SMD  *" $line]} { continue }
	    
	    set tok [concat $line]
	    set pos [lrange $tok 2 4]
	    set force [lrange $tok 5 7]

	    incr readCount
	    if {$readCount <= $skipSteps} {
		set pos0 $pos
		continue 
	    }
	    set line [makeVelocityForce grid $pos0 $pos $force [expr {$dt*$readCount}] $dt]
	    set pos0 $pos

	    if {[llength $lineList] >= $trajSteps} {
		# Find the node of this instance.
		# Watch for the periodic boundaries.
		set cen [wrap grid [meanPos [wrapToSelf grid $posList]]]
		set node [cylGeoNearestIndex cyl $cen]

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

proc distributeCorrCylBinary {gridVar cylVar inFileList outDir outName skipSteps trajSteps dt} {
    upvar $gridVar grid
    upvar $cylVar cyl
    set timestep 1.0

    foreach inFile $inFileList {
	# Open the files.
	set in [open $inFile r]
	set inData [open $inFile r]
	puts "Reading `$inFile...'"

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
	    set line [makeVelocityForce grid $pos0 $pos1 $force [expr {$dt*$readCount}] $timestep]
	    set inBinData [string range $inBinData $chunk end]

	    if {[llength $lineList] >= $trajSteps} {
		# Find the node of this instance.
		# Watch for the periodic boundaries.
		set cen [wrap grid [meanPos [wrapToSelf grid $posList]]]
		set node [cylGeoNearestIndex cyl $cen]

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

# Get the cylinder geometry.
cylGeoRead cyl $cylFile
puts "Read cylindrical grid: "
puts "s0 $cyl(s0) z0 $cyl(z0)"
puts "ds $cyl(ds) dz $cyl(dz)"
puts "ns $cyl(ns) nz $cyl(nz)"

# Get the grid geometry.
readDx grid $gridFile
puts "Read the grid with dimensions $grid(nx) $grid(ny) $grid(nz)."

# Delete the files in the output directory that would collide.
set collideList [glob -nocomplain $outDir/$outName.*.dat]
if {[llength $collideList] > 0} {
    foreach collide $collideList {
	#file delete $collide
	#puts "Deleted $collide"	
    }
}

# Sequester each instance (mini-trajectory) in a file for the corresponding grid node.
if {$isBinary} {
    set n [distributeCorrCylBinary grid cyl $inFileList $outDir $outName $skipSteps $instanceSteps $dt]
} else {
    set n [distributeCorrCyl grid cyl $inFileList $outDir $outName $skipSteps $instanceSteps $dt]
}
puts "Distributed $n instances."

exit
