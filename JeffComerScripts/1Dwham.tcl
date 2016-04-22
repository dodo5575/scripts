# whamPmf.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

proc computeSpring {springK r0 r} {
    set u 0.0
    foreach kx $springK x0 $r0 x $r {
	set u [expr $u + 0.5*$kx*($x-$x0)*($x-$x0)]
    }
    return $u
}

proc whamUsage {} {
    puts "Usage: wham <xVar> <winCen> <winSpring> <tol?>"
    puts "<xVar> is the name of a list of lists of coordinate data for each window."
    puts "<x> is indexed by (window step)."
    puts "<winCen> is a list of the positions of the harmonic springs for each window"
    puts "<winSpring> is a list of the spring constants for each window (in k_B*T/units(<x>)^2)"
    puts "<winSpring> can also be a scalar in the case where the spring constants are identical for all windows."
    puts "<tol> controls the number of WHAM iterations."
    puts "The result is considered sufficiently accurate when the maximum change in any F value (in k_B*T) over an iteration is less than <tol>."
    error ""
}

# Compute the PMF from a probability distribution.
proc generatePmf {binCen prob} {
    set n [llength $prob]

    # Find the maximum and minimum values.
    for {set b 0} {$b < $n} {incr b} {
	set p [lindex $prob $b]
			    
	if {$p > 0.0} { break }
    }

    if {$b == $n} {
	puts "Error! Probability distribution is zero."
	return
    }
    set pmfMin [expr -log($p)]
    set pmfMax [expr -log($p)]
    foreach p $prob {
	if {$p > 0.0} {
	    set u [expr -log($p)]
	    if {$u < $pmfMin} { set pmfMin $u }
	    if {$u > $pmfMax} { set pmfMax $u }
	}
    }
    
    # Calculate the PMF.
    # The PMF is associated with the bin centers.
    set pmf {}
    foreach c $binCen p $prob {
	if {$p > 0.0} {
	    lappend pmf [list $c [expr -log($p) - $pmfMin]]
	} else {
	    lappend pmf [list $c [expr $pmfMax - $pmfMin]]
	}
    }

    return $pmf
}

# Run the Weighted Histogram Analysis Method in one dimension.
# $xVar is the name of a list of lists of coordinate data for each window.
# $x is indexed by (window step).
# $winCen is list of window centers indexed by (window).
# $winSpring is the spring constant used in every window in units of k_B*T/units($x)^2.
# $binN is the number of bins along each direction.
# $tol is the convergence energy in k_B*T.
#
# See: B. Roux. The calculation of the potential of mean force using computer
# simulations. Computer Physics Communications 91, 275 (1995).
proc wham {xVar binN winCen winSpring outFile levelX0 levelX1 {tol 0.001} {type "linear"}} {
    upvar $xVar x
    set energyThreshold 1e20

    # Validate the input.
    if {[llength $winCen] != [llength $x]} {
	puts "Error:wham Inconsistent dimensions x([llength $x]) and winCen([llength $winCen])."
	error ""
    }

    # Make a list out of $winSpring if it is a scalar.
    if {[llength $winSpring] == 1} {
	set ws $winSpring
	set winSpring {}
	foreach w $winCen {
	    lappend winSpring $ws
	}
    } elseif {[llength $winSpring] != [llength $winCen]} {
	puts "Error:wham Inconsistent dimensions winSpring([llength $winSpring]) and winCen([llength $winCen])."
	error ""
    }

    set Pi [expr {4*atan(1)}]

    puts "Using Pi = $Pi"

    # Form the bins.
    foreach {binOrigin binSize} [makeBins0 x $binN] {break}
    puts "Bin origin: $binOrigin"
    puts "Bin size: $binSize"
    puts "Number of bins: $binN"

    # Determine the bin volume -- for bond lengths, it's a shell, not a line
    array set volume [list ]
    for {set i 0} {$i < $binN} {incr i} {
	if { [string equal $type "linear"] } {
	    set volume($i) $binSize
	} elseif { [string equal $type "shell"] } {
	    set r0 [expr {$binOrigin + $i*$binSize}]
	    set r1 [expr {$binOrigin + ($i+1)*$binSize}]
	    set volume($i) [expr {(4*$Pi / 3.0)*($r1*$r1*$r1-$r0*$r0*$r0)}]
	    puts "volume($i) is $volume($i)"
	} else {
	    print "There is no defined volume of type $type\nExiting...."
	    return 0;
	}
    }
    
    # Histogram the data to form the biased probability distributions for each window.
    # $winProbBias is indexed by (window bin)
    # $winSteps is indexed by (window)
    set outHist [open $outFile.hist w]
    set winProbBias {}
    set winSteps {}
    set winList {}
    set winHist {}
    set w 0
    foreach win $x c $winCen ws $winSpring {
	set histN [histogram win count $binOrigin $binSize $binN]
	#lappend winSteps [llength $win]
	lappend winSteps $histN
	puts "Histogrammed $histN of [llength $win] data points for window at $c."

	# Compute the average energy for the window.
	set energy 0.0
	foreach r $win { set energy [expr $energy + 0.5*$ws*($c-$r)*($c-$r)] }
	if {[llength $win] > 0} {
	    set energy [expr $energy/[llength $win]]
	    puts "Average energy is $energy kT."
	}

	# Decide whether to keep the window or not.
	if {$histN > 0 && [llength $win] > 0 && $energy < $energyThreshold} {
	    lappend winList $w
	    
	    set freq {}
	    set hist {}
	    for {set i 0} {$i < $binN} {incr i} {
		puts -nonewline $outHist "$count($i) "
		lappend hist $count($i)
		
		lappend freq [expr {double($count($i))/$histN/$volume($i)}]
	    }

	    lappend winHist $hist
	    lappend winProbBias $freq
	    puts $outHist ""
	}
	incr w
    }
    close $outHist
    
    # Make a list of bin centers.
    set binCen {}
    for {set i 0} {$i < $binN} {incr i} {
	lappend binCen [expr {$binOrigin + ($i+0.5)*$binSize}]
    }
    puts "Bin centers: [llength $binCen]"

    
    # Write total histogram.
    set outHist [open $outFile.histTotal w]
    for {set b 0} {$b < $binN} {incr b} {
	set sum 0.0
	foreach w $winList hist $winHist {
	    # Get the biased prob. distrib. at this bin center for this window.
	    set sum [expr $sum + [lindex $hist $b]]
	}
	puts $outHist "[lindex $binCen $b] $sum"
    }
    close $outHist

    # Make an initial guess for the f constants. 
    set winF {}
    foreach w $winList {
	lappend winF 0.0
    }
    set outF [open $outFile.fConst w]

    # Compute the numerator of Roux Eq. 8 beforehand.
    set binNumer {}
    for {set b 0} {$b < $binN} {incr b} {
	set numer 0.0; # numerator of Roux Eq. 8.

	foreach w $winList probBias $winProbBias {
	    set n [lindex $winSteps $w]

	    # Get the biased prob. distrib. at this bin center for this window.
	    set p [lindex $probBias $b]

	    # Add to the numerator sum in Roux Eq. 8.
	    set numer [expr {$numer + $n*$p}]
	}
	lappend binNumer $numer
    }

    # Iterate the WHAM method.
    set maxChange 1
    set iter 0
    while {$maxChange > $tol} { 
	# Roux Eq. 8
	#
	# Estimate the unbiased probability distribution.
	# Loop through the bins.
	set prob {}; # Estimate of unbiased prob. distrib. Indexed by (bin).
	foreach bc $binCen numer $binNumer {
	    # Compute the probability distribution at this bin center. 
	    set denom 0.0; # denominator of Roux Eq. 8.
	    
	    # Loop through the windows to do the sums in Roux Eq. 8.
	    foreach w $winList f $winF {
		set n [lindex $winSteps $w]
		set wc [lindex $winCen $w]
		set ws [lindex $winSpring $w]

	  	# Compute the bias energy at this bin center for this window.
		set dr [expr {$bc - $wc}]
		set uBias [expr {0.5*$ws*$dr*$dr}]

		# Add to the denomenator sum in Roux Eq. 8.
		set denom [expr {$denom + $n*exp($f-$uBias)}]
	    }

	    # Compute the unbiased prob. distrib. for this bin.
	    lappend prob [expr {$numer/$denom}]
	}

	# Roux Eq. 9
	#
	# Integrate to obtain the f constants for each window.
	set winF0 $winF
	set winF {}
	foreach w $winList {
	    set c [lindex $winCen $w]
	    set ws [lindex $winSpring $w]
	    
	    # Numerically integrate over the reaction coordinate (over the bins).
	    set sum 0.0
	    for {set b 0} {$b < $binN} {incr b} {
		set dr [expr {[lindex $binCen $b]-$c}]
		set uBias [expr {0.5*$ws*$dr*$dr}]
	        set sum [expr {$sum + $volume($b)*exp(-$uBias)*[lindex $prob $b]}]
	    }
	    
	    set f [expr {-log($sum)}]
	    lappend winF $f
	}

	set winFTemp {}
	foreach f $winF {
	    lappend winFTemp [expr {$f - [lindex $winF end]}]
	}
	#set winF $winFTemp

	# Check for convergence.
	# Compute the differences between consecutive f values.
	# Find the largest change in these differences.
	set w 1
	set df0 [expr {[lindex $winF0 $w]-[lindex $winF0 [expr $w-1]]}]
	set df [expr {[lindex $winF $w]-[lindex $winF [expr $w-1]]}]
	set maxChange [expr {abs($df-$df0)}]
	for {set w 2} {$w < [llength $winF]} {incr w} {
	    set df0 [expr {[lindex $winF0 $w]-[lindex $winF0 [expr $w-1]]}]
	    set df [expr {[lindex $winF $w]-[lindex $winF [expr $w-1]]}]
	    if {[expr {abs($df-$df0)}] > $maxChange} {
		set maxChange [expr {abs($df-$df0)}]
	    }
	}
	puts "WHAM iteration $iter: $maxChange"
	incr iter
	puts $outF $winF
    }
    close $outF

    set pmf [generatePmf $binCen $prob]

    # --- Find the level value, and the error, stddev --- #
    set levelSum 0.0
    set levelCount 0
    foreach pair $pmf {
	set z [lindex $pair 0]
	set u [lindex $pair 1]
	
	if {$z >= $levelX0 && $z < $levelX1} {
	    set levelSum [expr {$levelSum + $u}]
	    incr levelCount
	}
    }

    if {$levelCount > 2} {
	set levelMean [expr {$levelSum/$levelCount}]
	set levelSum 0.0
	set levelCount 0
	foreach pair $pmf {
	    set z [lindex $pair 0]
	    set u [lindex $pair 1]
	    
	    if {$z >= $levelX0 && $z < $levelX1} {
		set levelSum [expr {$levelSum + ($levelMean-$u)*($levelMean-$u)}]
		incr levelCount
	    }
	}
	set levelError [expr {1.96*sqrt($levelSum)/$levelCount}]
	set levelStdDev [expr {sqrt($levelSum/($levelCount-1))}]
	
	puts "Zero level: $levelMean +- $levelError / $levelStdDev"
    } else {
	puts "No data between $levelX0 and $levelX1."
	set levelMean 0.0
    }

    # --- --- --- --- --- --- --- --- --- --- --- --- --- #

    set out [open $outFile.dat w]
    foreach pair $pmf {
	puts $out "[lindex $pair 0] [expr [lindex $pair 1] - $levelMean] $levelStdDev"
    }
    close $out
    
    # We are finished with the WHAM iterations.
    return $pmf
}


# Run the Weighted Histogram Analysis Method in three dimensions.
# $rVar is the name of a list of lists of coordinate data for each window.
# $r is indexed by (window step).
# $winCen is list of window centers (vector) indexed by (window).
# $winSpring is the vector spring constant used in every window in units of k_B*T/units($x)^2.
# $bins is a vector with the number of bins along each direction.
# $tol is the convergence energy k_B*T.
#
# See: B. Roux. The calculation of the potential of mean force using computer
# simulations. Computer Physics Communications 91, 275 (1995).
proc wham3d {rVar winCen winSpring bins {tol 0.001}} {
    upvar $rVar r

    # Validate the input.
    if {[llength $winCen] != [llength $r]} {
	error "Error:wham3d Inconsistent dimensions r([llength $r]) and winCen([llength $winCen])."
	exit
    }
    
    if {[llength [lindex $winCen 0]] != 3} {
	error "Error:wham3d window centers \$winCen must be 3-vectors."
    }

    if {[llength [lindex $winSpring 0]] != 3 && [llength $winSpring] != 3} {
	error "Error:wham3d window spring constants \$winSpring  must be 3-vectors."
    }

    # Make a list out of $winSpring if it is only one vector.
    if {[llength $winSpring] == 3} {
	set ws $winSpring
	set winSpring {}
	foreach w $winCen {
	    lappend winSpring $ws
	}
    } elseif {[llength $winSpring] != [llength $winCen]} {
	puts "Error:wham3d Inconsistent dimensions winSpring([llength $winSpring]) and winCen([llength $winCen])."
	error ""
    }
    puts "\nComputing the pmf by the WHAM method with [llength $winCen] windows."
    puts "Based on B. Roux. Computer Physics Communications 91, 275 (1995)."

    # Form the bins.
    set binNx [lindex $bins 0]
    set binNy [lindex $bins 1]
    set binNz [lindex $bins 2]
    set binN [expr $binNx*$binNy*$binNz]
    set minMax [getMinMax3d2 $r]
    set binOrigin [lindex $minMax 0]
    set binSize {}
    set binVol 1.0
    foreach x0 [lindex $minMax 0] x1 [lindex $minMax 1] n $bins {
	set l [expr ($x1-$x0)/$n]
	lappend binSize $l
	set binVol [expr $binVol*$l]
    }
    puts "\nUsing $binN bins"
    puts "Bin origin: $binOrigin"
    puts "Bin size: $binSize"
    puts "Bin volume: $binVol"
    
    # Make a list of bin centers.
    set binCen {}
    for {set iz 0} {$iz < $binNz} {incr iz} {
	for {set iy 0} {$iy < $binNy} {incr iy} {
	    for {set ix 0} {$ix < $binNx} {incr ix} {
		set x [expr [lindex $binOrigin 0] + ($ix+0.5)*[lindex $binSize 0]]
		set y [expr [lindex $binOrigin 1] + ($iy+0.5)*[lindex $binSize 1]]
		set z [expr [lindex $binOrigin 2] + ($iz+0.5)*[lindex $binSize 2]]
		
		lappend binCen [list $x $y $z]
	    }
	}
    }

    # Histogram the data to form the biased probability distributions for each window.
    # $winProbBias is indexed by (window bin)
    # $winSteps is indexed by (window)
    puts "Histogramming data..."
    set winProbBias {}
    set winSteps {}
    foreach win $r wc $winCen {
	set histN [histogram3d win count $binOrigin $binSize $bins]
	lappend winSteps $histN
	 puts "Histogrammed $histN data points data for the window with center $wc."

	set freq {}
	for {set i 0} {$i < $binN} {incr i} {
	    lappend freq [expr double($count($i))/[llength $win]/$binVol]
	}

	lappend winProbBias $freq
    }

    # Compute the numerator of Roux Eq. 8 beforehand.
    puts -nonewline "\nPrecomputing the numerator of Roux Eq. 8..."
    set binNumer {}
    for {set b 0} {$b < $binN} {incr b} {
	set numer 0.0; # numerator of Roux Eq. 8.

	foreach probBias $winProbBias n $winSteps {
	    # Get the biased prob. distrib. at this bin center for this window.
	    set p [lindex $probBias $b]

	    # Add to the numerator sum in Roux Eq. 8.
	    set numer [expr $numer + $n*$p]
	}
	lappend binNumer $numer
    }
    puts "Done."

    # Compute exp(-uBias) beforehand.
    puts -nonewline "Precomputing exp(-uBias) for Roux Eq. 8 and 9..."
    set bi 0
    foreach bc $binCen {
	set wi 0
	foreach wc $winCen ws $winSpring {
	    set uBias [computeSpring $ws $wc $bc]
	    set biasFactor($bi,$wi) [expr exp(-$uBias)]
	    incr wi
	}
	incr bi
    }
    puts "Done."

    # Make an initial guess for the f constants. 
    set winF {}
    foreach c $winCen {
	lappend winF 0.0
    }

    # Iterate the WHAM method.
    puts "\nBeginning WHAM iterations."
    set maxChange 1
    set iter 0
    while {$maxChange > $tol} { 
	# Roux Eq. 8
	#
	# Estimate the unbiased probability distribution.
	# Loop through the bins.
	set prob {}; # Estimate of unbiased prob. distrib. Indexed by (bin).
	set bi 0
	foreach bc $binCen numer $binNumer {
	    # Compute the probability distribution at this bin center. 
	    set denom 0.0; # denominator of Roux Eq. 8.
	    
	    # Loop through the windows to do the sums in Roux Eq. 8.
	    set wi 0
	    foreach n $winSteps f $winF wc $winCen {
		# Compute the bias energy at this bin center for this window.
		# Add to the denomenator sum in Roux Eq. 8.
		
		set denom [expr $denom + $n*exp($f)*$biasFactor($bi,$wi)]
		incr wi
	    }

	    # Compute the unbiased prob. distrib. for this bin.
	    lappend prob [expr $numer/$denom]
	    incr bi
	}
	
	# Roux Eq. 9
	#
	# Integrate to obtain the f constants for each window.
	set winF0 $winF
	set winF {}
	set wi 0
	foreach wc $winCen {
	    # Numerically integrate over the reaction coordinate (over the bins).
	    set sum 0.0
	    for {set bi 0} {$bi < $binN} {incr bi} {
	        set sum [expr $sum + $binVol*$biasFactor($bi,$wi)*[lindex $prob $bi]]
	    }
	    
	    set f [expr -log($sum)]
	    lappend winF $f
	    incr wi
	}

	# Check for convergence.
	# Compute the differences between consecutive f values.
	# Find the largest change in these differences.
	set w 1
	set df0 [expr [lindex $winF0 $w]-[lindex $winF0 [expr $w-1]]]
	set df [expr [lindex $winF $w]-[lindex $winF [expr $w-1]]]
	set maxChange [expr abs($df-$df0)]
	for {set w 2} {$w < [llength $winF]} {incr w} {
	    set df0 [expr [lindex $winF0 $w]-[lindex $winF0 [expr $w-1]]]
	    set df [expr [lindex $winF $w]-[lindex $winF [expr $w-1]]]
	    if {[expr abs($df-$df0)] > $maxChange} {
		set maxChange [expr abs($df-$df0)]
	    }
	}
	puts "WHAM iteration $iter: $maxChange"
	incr iter
    }

    # We are finished with the WHAM iterations.
    # Find the minimum and maximum of the pmf.
    set pmfMin 1e100
    set pmfMax -1e100
    foreach p $prob {
	if {$p > 0.0} {
	    set u [expr -log($p)]
	   if {$u < $pmfMin} {set pmfMin $u}
	   if {$u > $pmfMax} {set pmfMax $u}
	}
    }

    # Form the pmf list with the minimum zero.
    # Replace infinite potentials with the maximum value.
    set pmf {}
    foreach p $prob {
	if {$p > 0.0} {
	    set u [expr -log($p)]
	    lappend pmf [expr $u - $pmfMin]
	} else {
	    lappend pmf [expr $pmfMax-$pmfMin]
	}
    }

    return [list $binCen $pmf]
}

# Return the minimum and maximum values in a list.
proc getMinMax {valList} {
    set minVal [lindex $valList 0] 
    set maxVal [lindex $valList 0] 
    foreach x $valList {
	if {$x < $minVal} {set minVal $x}
	if {$x > $maxVal} {set maxVal $x}
    }
    return [list $minVal $maxVal]
}

# Return the minimum and maximum values in a list of lists.
proc getMinMax2 {valList} {
    set minVal [lindex $valList 0 0] 
    set maxVal [lindex $valList 0 0] 
    foreach v $valList {
	foreach x $v {
	    if {$x < $minVal} {set minVal $x}
	    if {$x > $maxVal} {set maxVal $x}
	}
    }
    return [list $minVal $maxVal]
}

# Return the minimum and maximum values in a list of lists.
proc getMinMax3d2 {valList} { 
    set x0 [lindex $valList 0 0 0]
    set x1 [lindex $valList 0 0 0] 
    set y0 [lindex $valList 0 0 1]
    set y1 [lindex $valList 0 0 1] 
    set z0 [lindex $valList 0 0 2]
    set z1 [lindex $valList 0 0 2] 

    foreach v $valList {
	foreach r $v {
	    foreach {x y z} $r {break}
	    if {$x < $x0} {set x0 $x}
	    if {$x > $x1} {set x1 $x}
	    if {$y < $y0} {set y0 $y}
	    if {$y > $y1} {set y1 $y}
	    if {$z < $z0} {set z0 $z}
	    if {$z > $z1} {set z1 $z}
	}
    }
    return [list [list $x0 $y0 $z0] [list $x1 $y1 $z1]]
}

# Using the range of the list with the name $xVar, make $nx
# bins.
# Returns the starting location and the size of the $nx bins.
proc makeBins0 {xVar nx} {
    upvar $xVar x

    set sx [getMinMax2 $x]
    foreach {x0 x1} $sx { break }
    set dx [expr ($x1-$x0)/$nx]

    return [list $x0 $dx]
}

proc makeBins {xVar nx factor0 factor1} {
    upvar $xVar x

    set sx [getMinMax2 $x]
    set lx [expr [lindex $sx 1] - [lindex $sx 0]]
    set x0 [expr [lindex $sx 0] + $factor0*$lx]
    set x1 [expr [lindex $sx 1] - $factor1*$lx]
    set dx [expr ($x1-$x0)/$nx]

    return [list $x0 $dx]
}

proc makeBin {min max n factor0 factor1} {
    set lx [expr $max - $min]
    set x0 [expr $min + $factor0*$lx]
    set x1 [expr $max - $factor1*$lx]
    set dx [expr ($x1-$x0)/$n]

    return [list $x0 $dx]
}

# Histogram some data.
# $rxVar is the name of list containing the data.
# $countVar will be an array with counts in each bin.
# $x0 is the origin of the bins.
# $dx is the size of the bins.
# $n is the number of bins.
# Returns the starting position and the size of the $n bins.

proc histogram {rxVar countVar x0 dx n} {
    upvar $rxVar rx
    upvar $countVar count
    set samples 0

    # Zero the bin counts and form the bin ce
    for {set j 0} {$j < $n} {incr j} {
	set count($j) 0
    }    

    # Add the bin counts.
    foreach x $rx {
	set ix [expr int(floor(($x - $x0)/$dx))]
	if {$ix >= 0 && $ix < $n} {
	    incr count($ix)
	    incr samples
	}
    }
    return $samples
}


proc histogram3d {rVar countVar origin delta nodes} {
    upvar $countVar count
    upvar $rVar r
    set samples 0
  
    foreach {nx ny nz} $nodes {break}
    foreach {x0 y0 z0} $origin {break}
    foreach {dx dy dz} $delta {break}
    set n [expr $nx*$ny*$nz]
    set nxny [expr $nx*$ny]

    # Zero the bin counts.
    for {set j 0} {$j < $n} {incr j} {
	set count($j) 0
    }
    
    # Add the bin counts.
    foreach p $r {
	foreach {x y z} $p {break}
	set ix [expr int(floor(($x - $x0)/$dx))]
	set iy [expr int(floor(($y - $y0)/$dy))]
	set iz [expr int(floor(($z - $z0)/$dz))]

	if {$ix < 0} {continue}
	if {$ix >= $nx} {continue}
	if {$iy < 0} {continue}
	if {$iy >= $ny} {continue}
	if {$iz < 0} {continue}
	if {$iz >= $nz} {continue}

	set j [expr $ix + $nx*$iy + $nxny*$iz]
	incr count($j)
	incr samples
    }
    return $samples
}

# Just read a space delimited data file.
proc readData {fileName} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] < 2} {continue}

	set tok [concat $line]
	set t [lindex $tok 0]
	lappend r $tok
    }

    close $in
    return $r
}

# Read the position data, ignoring data younger than $cutTime.
proc readPositions {fileName cutTime} {
    set in [open $fileName r]
    
    set r {}
    while {[gets $in line] >= 0} {
	if {[string match "#*" $line]} {continue}
	if {[string length $line] <= 1} {continue}

	set tok [concat $line]
	set t [lindex $tok 0]
	if {$t > $cutTime} {
	    lappend r [lrange $tok 1 3]
	}
    }

    close $in
    return $r
}

proc reverseIndex {lVar nodes} {
    upvar $lVar l
    set ret {}
    
    foreach {nx ny nz} $nodes {break}
    for {set ix 0} {$ix < $nx} {incr ix} {
	for {set iy 0} {$iy < $ny} {incr iy} {
	    for {set iz 0} {$iz < $nz} {incr iz} {
		set j [expr $ix + $nx*$iy + $nx*$ny*$iz]

		lappend ret [lindex $l $j]
	    }
	}
    }
    return $ret
}

# The potential value at each point is in gridVar(data),
# the number of grid points along each
# cell axis is in gridVar(nx,ny,nz), the lattice vectors are in
# gridVar(delta) (3x3 matrix), and the cell origin is in
# gridVar(origin) (3-vector).
# Data order is z fast, y medium, and x slow.
# gridVar(size) = nx*ny*nz.
# gridVar(deltaInv) is the inverse of the delta matrix.
# Writes a dx file given a grid structure.
proc writeDx {gridVar fileName} {
    upvar $gridVar grid

    foreach {ox oy oz} $grid(origin) {break}
    set delta [join $grid(delta)]
    foreach {exx eyx ezx exy eyy ezy exz eyz ezz} $delta {break}

    set out [open $fileName w]
    # Write the headers.
    puts $out "\# NAMD gridforce grid"
    # Write the grid attributes.
    puts $out "object 1 class gridpositions counts $grid(nx) $grid(ny) $grid(nz)"
    puts $out "origin $ox $oy $oz"
    puts $out "delta $exx $exy $exz"
    puts $out "delta $eyx $eyy $eyz"
    puts $out "delta $ezx $ezy $ezz"
    puts $out "object 2 class gridconnections counts $grid(nx) $grid(ny) $grid(nz)"
    puts $out "object 3 class array type double rank 0 items $grid(size) data follows"
    
    # Write the data.
    set penultima [expr 3*($grid(size)/3)]
    set mod [expr $grid(size)-$penultima]
    set i 0
    foreach {v0 v1 v2} $grid(data) {
	puts $out "$v0 $v1 $v2"
	
	incr i 3
	# Treat the last line specially.
	if {$i >= $penultima} {break}
    }

    if {$mod == 1} {
	puts $out "[lindex $grid(data) end]"
    } elseif {$mod == 2} {
	puts $out "[lindex $grid(data) end-1] [lindex $grid(data) end]"
    }

    close $out
}

# Writes the dx file header given an output stream.
proc writeDxHeader {gridVar out} {
    upvar $gridVar grid

    foreach {ox oy oz} $grid(origin) {break}
    set delta [join $grid(delta)]
    foreach {exx eyx ezx exy eyy ezy exz eyz ezz} $delta {break}

    # Write the headers.
    puts $out "\# NAMD gridforce grid"
    # Write the grid attributes.
    puts $out "object 1 class gridpositions counts $grid(nx) $grid(ny) $grid(nz)"
    puts $out "origin $ox $oy $oz"
    puts $out "delta $exx $exy $exz"
    puts $out "delta $eyx $eyy $eyz"
    puts $out "delta $ezx $ezy $ezz"
    puts $out "object 2 class gridconnections counts $grid(nx) $grid(ny) $grid(nz)"
    puts $out "object 3 class array type double rank 0 items $grid(size) data follows"
}
