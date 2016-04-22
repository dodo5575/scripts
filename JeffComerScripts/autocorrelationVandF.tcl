#!/usr/bin/tclsh
# Compute 
# Author: Jeff Comer <jcomer2@illinois.edu>

source vector.tcl
source gridForce.tcl

if {$argc != 3} {
    puts "Usage: ./correlationFunctions.tcl gridFile inputFile0 outName"
    exit
}

set gridFile [lindex $argv 0]
set inFileList [list [lindex $argv 1]]
set outName [lindex $argv 2]
set dt 1.0
set kT 0.5862292
set doVelocity 1
set doForce 1
set minSamples 6

proc meanTrajList {listList} {
    set meanList {}
    set n [llength [lindex $listList 0]]
    
    for {set i 0} {$i < $n} {incr i} {
	set sum 0.0
	set count [llength $listList]
	foreach vec $listList {
	    set sum [expr {$sum + [lindex $vec $i]}]
	}
	lappend meanList [expr $sum/$count]
    }
    
    return $meanList
}

proc meanTrajListVec {listList} {
    set meanList {}
    set n [llength [lindex $listList 0]]
    
    for {set i 0} {$i < $n} {incr i} {
	set sum [vecZero]
	set count [llength $listList]
	foreach vec $listList {
	    set sum [vecAdd $sum [lindex $vec $i]]
	}
	lappend meanList [vecScale [expr {1.0/$count}] $sum]
    }
    
    return $meanList
}

proc addTrajListVec {listList vectorList scale} {
    set retTrajList {}
    foreach lis $listList {
	set retList {}
	foreach vector $lis vector0 $vectorList {
	    lappend retList [vecAdd $vector [vecScale $scale $vector0]]
	}
	lappend retTrajList $retList
    }
    return $retTrajList
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

proc loadInstances {instancePosVar instanceForceVar inFileList} {
    upvar $instancePosVar instancePos
    upvar $instanceForceVar instanceForce

    foreach inFile $inFileList {
	# Open the files.
	set in [open $inFile r]
	puts "Read $inFile."

	set posList {}
	set forceList {}

	set count 0
	while {[gets $in line] >= 0} {
	    if {[string length $line] <= 1} { continue }
	    if {[string match "#*" $line]} { continue }
	    
	    # We've reached the end of an instance. Store it.
	    if {[string match "END*" $line]} {
		# Store the position and force lists.
		lappend instancePos $posList
		lappend instanceForce $forceList

		set posList {}
		set forceList {}
		incr count
		continue
	    }

	    set item [concat $line]
	    if {[llength $item] < 7} { 
		puts "Only [llength $item] of 7 columns found."
		continue 
	    }
	    
	    # Build up the lists.
	    lappend posList [lrange $line 1 3]
	    lappend forceList [lrange $line 4 6]
	}
	close $in
    }
    return
}

proc differentiate {posList dt} {
    set velList {}
    set n [llength $posList]

    for {set i 1} {$i < $n} {incr i} {
	set v [vecScale [expr {1.0/$dt}] [vecSub [lindex $posList $i] [lindex $posList [expr $i-1]]]]
	lappend velList $v
    }
    return $velList
}

proc correlationVec {valList} {
    set corr {}
    set v0 [lindex $valList 0]
    foreach v $valList {
	lappend corr [vecDot $v0 $v]
    }
    return $corr
}

# Get the grid geometry.
readDx grid $gridFile
puts "Read the grid with dimensions $grid(nx) $grid(ny) $grid(nz)."


# Load the instances (mini-trajectories).
set instPosList {}
set instForceList {}
loadInstances instPosList instForceList $inFileList
puts "Loaded [llength $instPosList] instances."
copyGridDim grid meanForceGrid
copyGridDim grid diffuseVelGrid
copyGridDim grid diffuseForceGrid
copyGridDim grid instCountGrid

# Compute the center for each trajectory.
# Assign it to a place in the grid.
set occupiedGridPoints {}
set i 0
foreach pos $instPosList {
    set cen [vecZero]
    # Compute the mean position.
    foreach r $pos {
	set cen [vecAdd $cen $r]
    }
    set cen [vecScale [expr {1.0/[llength $pos]}] $cen]

    # Use the mean position to set the associated grid point.
    set ind [nearestIndex grid $cen]

    # Make an array entry for each occupied grid point.
    if {[info exists gridInst($ind)]} {
	# The array entry already exists. Append its index.
	lappend gridInst($ind) $i
    } else {
	# This is a new array entry.
	set gridInst($ind) [list $i]
	lappend occupiedGridPoints $ind
    }
    lset instCountGrid(data) $ind [expr {[lindex $instCountGrid(data) $ind] + 1}]
    incr i
}
set instOutFile ${outName}_instances.dx
writeDx instCountGrid $instOutFile
puts "Wrote the number of instance counts applied to each grid point to $instOutFile."

# Compute the velocity correlation function for each grid point.
if {$doVelocity} {
    set corrVelAll {}
    foreach occIndex $occupiedGridPoints {
	set instList $gridInst($occIndex)
	
	# Compute the velocity correlation functions.
	set corrTrajList {}
	foreach inst $instList {
	    # Get the position list from the instance list.
	    set posTraj [lindex $instPosList $inst]
	    set velTraj [differentiate $posTraj $dt]
	    set correlationFunc [correlationVec $velTraj]
	    lappend corrTrajList $correlationFunc
	    lappend corrVelAll $correlationFunc
	}

	# Compute the mean velocity correlation function at this grid point.
	set meanCorrVel [meanTrajList $corrTrajList]
	# Write it.
	#writeCorrFunction ${outName}${occIndex}_vel.dat $meanCorrVel $dt

	# Integrate the correlation function to obtain the diffusivity.
	set sum 0.0
	foreach corr $meanCorrVel {
	    set sum [expr $sum + $corr*$dt]
	}
	# Scale by the number of dimensions.
	set diffuse [expr {$sum/3.0}]
	# Convert to A^2/ns. Was in A^2/fs.
	set diffuse [expr {$diffuse*1e6}]

	# Set this point in the grid.
	# Leave it at zero if there aren't enough samples.
	if {[lindex $instCountGrid(data) $occIndex] >= $minSamples} {
	    lset diffuseVelGrid(data) $occIndex $diffuse
	}
    }

    set meanCorrVelAll [meanTrajList $corrVelAll]
    set outFile ${outName}_velCorr.dat
    writeCorrFunction $outFile $meanCorrVelAll $dt
    puts "Wrote average force correlation function for the entire system to $outFile."

    set velOutFile ${outName}_vel.dx
    writeDx diffuseVelGrid $velOutFile
    puts "Wrote the diffusion map derived from the velocity autocorrelation function to $velOutFile."
}

# Compute the force correlation function for each grid point.
if {$doForce} {
    set corrForceAll {}
    foreach occIndex $occupiedGridPoints {
	set instList $gridInst($occIndex)
	
	# Collect the forces.
	set forceTrajList {}
	foreach inst $instList {
	    lappend forceTrajList [lindex $instForceList $inst]
	}

	# Get the mean force over the trajectories <F(r,t)>.
	set meanForce [meanTrajListVec $forceTrajList]
	lset meanForceGrid(data) $occIndex $meanForce

	# Compute the deviation of the force F(r,t) - <F(r,t)>.
	set deltaForceTrajList [addTrajListVec $forceTrajList $meanForce -1.0]
	
	# Compute the force correlation functions.
	set corrTrajList {}
	foreach deltaForceTraj $deltaForceTrajList {
	    set correlationFunc [correlationVec $deltaForceTraj]
	    lappend corrTrajList $correlationFunc
	    lappend corrForceAll $correlationFunc
	}

	# Compute the mean force correlation function at this grid point.
	set meanCorrForce [meanTrajList $corrTrajList]
	# Write it.
	#writeCorrFunction ${outName}${occIndex}_force.dat $meanCorrForce $dt

	# Integrate the correlation functions to obtain the diffusivity.
	set sum 0.0
	foreach corr $meanCorrForce {
	    set sum [expr $sum + $corr*$dt]
	}
	# Compute the diffusivity.
	if {$sum == 0.0} { continue }
	set diffuse [expr {$kT*$kT/$sum}]
	# Convert to A^2/ns. Was in A^2/fs.
	set diffuse [expr {$diffuse*1e6}]

	# Set this point in the grid.
	if {[lindex $instCountGrid(data) $occIndex] >= $minSamples} {
	    lset diffuseForceGrid(data) $occIndex $diffuse
	}
    }
    
    set meanCorrForceAll [meanTrajList $corrForceAll]
    set outFile ${outName}_forceCorr.dat
    writeCorrFunction $outFile $meanCorrForceAll $dt
    puts "Wrote average force correlation function for the entire system to $outFile."

    set meanForceOutFile ${outName}_meanForce.dx
    writeDx meanForceGrid $meanForceOutFile
    puts "Wrote the mean force at each grid point to $meanForceOutFile."

    set forceOutFile ${outName}_force.dx
    writeDx diffuseForceGrid $forceOutFile
    puts "Wrote the diffusion map derived from the force autocorrelation function to $forceOutFile."
}

puts "Complete."
exit
