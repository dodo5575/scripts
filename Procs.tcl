# proc append { list string } {
#     set tmp [join $list "$string "]
#     return "$tmp$string"
# }

source ~/scripts/vector.tcl

proc copy_cell_basis { molid_src molid_dst } {
    foreach prop {a b c alpha beta gamma} {
	molinfo $molid_dst set $prop [molinfo $molid_src get $prop]
    }
}

proc save_viewpoint { molid filename } {
    set ch [open $filename w]
    set matrices [molinfo $molid get {center_matrix rotate_matrix scale_matrix global_matrix}]
    puts $ch $matrices
    close $ch
}

proc restore_viewpoint { molids filename } {
    set ch [open $filename r]
    set matrices [gets $ch]
    foreach molid $molids {
	molinfo $molid set {center_matrix rotate_matrix scale_matrix global_matrix} $matrices
    }
    close $ch
}

proc selvel { sel } {
    # returns velocity of selection as determined by backward difference in time, i.e. looks at now-(now-1)
    # returns "0 0 0" if frame of $sel is already 0
    # return value is in Ang/frame
    set frame_init [$sel frame]
    if { $frame_init eq "now" } {
	set frame [molinfo top get frame]
    } else {
	set frame $frame_init
    }
    if { $frame == 0 } {
	set returnval {}
	for { set i 0 } { $i < [$sel num] } { incr i } {
	    lappend returnval {0 0 0}
	}
    } else {
	set r [$sel get {x y z}]
	$sel frame [expr {$frame - 1}]
	set r0 [$sel get {x y z}]
	set returnval [vecsub_list $r $r0]
    }
    $sel frame $frame_init
    return $returnval
}

proc mtime { file reffile deltat } {
    set mtime [expr [file mtime $reffile] - $deltat]
    file mtime $file $mtime
}

proc seldist { sel1 sel2 } {
    # measures distance between first atoms of sel1 and sel2 -- makes most sense for 1-atom selections
    set vec1 [lindex [$sel1 get {x y z}] 0]
    set vec2 [lindex [$sel2 get {x y z}] 0]
    return [veclength [vecsub $vec1 $vec2]]
}

proc fixions { {molid top} } {
    set ions	[list	LIT	SOD	POT	RUB	CES	CLA]
    set radii	[list	1.025	1.369	1.705	1.813	1.976	2.513]
    foreach ion $ions radius $radii {
	set sel [atomselect $molid "resname $ion"]
	$sel set radius $radius
	$sel delete
    }
}

proc wrap { sel args } {
    if { [llength $args] } {
	pbc wrap -sel $sel -compound residue -center origin $args
    } else {
	pbc wrap -sel $sel -compound residue -center origin -all
    }
}

proc wrapcom { sel sel2 args } {
    if { [llength $args] } {
	pbc wrap -sel $sel -compound residue -center com -centersel $sel2 $args
    } else {
	pbc wrap -sel $sel -compound residue -center com -centersel $sel2 -all
    }
}

proc com_align { {molid top} seltext } {
    set sel [atomselect $molid $seltext]
    set refsel [atomselect $molid $seltext frame 0]
    set refcom [measure center $refsel weight mass]
    set nframes [molinfo $molid get numframes]
    for { set frame 0 } { $frame < $nframes } { incr frame } {
	$sel frame $frame
	$sel moveby [vecsub $refcom [measure center $sel weight mass]]
    }
    $sel delete
    $refsel delete
}
	

proc water_density { {molid top} minmax } {
    lassign $minmax min max
    lassign $min xmin ymin zmin
    lassign $max xmax ymax zmax
    
    set vol [expr {($xmax-$xmin) * ($ymax-$ymin) * ($zmax-$zmin) / 1000.0}]
    set sel [atomselect $molid "name OH2 and x > $xmin and x < $xmax and y > $ymin and y < $ymax and z > $zmin and z < $zmax"]
    set nframes [molinfo $molid get numframes]
    
    #putvar xmin xmax ymin ymax zmin zmax
    
    set densities {}
    for { set frame 0 } { $frame < $nframes } { incr frame } {
	$sel frame $frame
	$sel update
	set num [$sel num]
	lappend densities [expr {$num/$vol}]
	puts [expr {$num/$vol}]
    }
    set density [lmean $densities]
    set density_stddev [lstddev $densities]
    puts "<density> = $density"
    puts "stddev(density) = $density_stddev"
}

proc putvar { args } {
    foreach var $args {
	upvar $var val
	puts "$var = $val"
    }
}

proc fputvar { ch args } {
    foreach var $args {
	upvar $var val
	puts $ch "# $var: $val"
    }
}

proc delete_all_sels { } {
    foreach sel [atomselect list] {
	#puts $sel
	uplevel $sel delete
    }
}

proc measurefit_zrot { sel1 sel2 } {
    # Like 'measure fit', but returns a matrix that only rotates about
    # the z axis. This is achieved by assuming 'measure fit' itself
    # returns a matrix for which this is nearly true already, then
    # computing the corresponding angles for the upper left 2x2
    # matrix, averaging these, and resetting the 3x3 matrix.
    
    # IMPORTANT NOTE: This works a lot better if the origin is
    # somewhere nearby, preferably at the center. The ignored
    # rotations are about the origin, so the further the origin, the
    # greater the error introduced. So the recommended procedure is to
    # move the two systems, find the matrix, apply it, then move back.
    
    set mat [measure fit $sel1 $sel2]
    
    set theta_xx [expr acos([lindex $mat 0 0])]
    set theta_xy [expr asin(-1 * [lindex $mat 0 1])]
    set theta_yx [expr asin(1 * [lindex $mat 1 0])]
    set theta_yy [expr acos([lindex $mat 1 1])]
    
    if { $theta_xy < 0 } {
	set theta_xx [expr -1 * $theta_xx]
	set theta_yy [expr -1 * $theta_yy]
    }
    
    set theta [expr ($theta_xx + $theta_xy + $theta_yx + $theta_yy)/4.0]
    puts "thetas = $theta_xx $theta_xy $theta_yx $theta_yy -> $theta"
    
    lset mat 0 0 [expr cos($theta)]
    lset mat 0 1 [expr -sin($theta)]
    lset mat 0 2 0.0
    lset mat 1 0 [expr sin($theta)]
    lset mat 1 1 [expr cos($theta)]
    lset mat 1 2 0.0
    lset mat 2 0 0.0
    lset mat 2 1 0.0
    lset mat 2 2 1.0
}

proc charge { sel } {
    measure sumweights $sel weight charge
}

proc basename { filename } {
    lindex [file split $filename] end
}

proc ionmolar { mol {resname CLA} } {
    # Finds molar concentration of ions
    # assumes there are as many cations as anions
    set selWat [atomselect $mol "name OH2"]
    set nWat [$selWat num]
    $selWat delete
    
    set selCLA [atomselect $mol "resname $resname"]
    set nCLA [$selCLA num]
    $selCLA delete
    
    puts "N(water) = $nWat"
    puts "N($resname) = $nCLA"
    
    expr 55.6 * $nCLA / ($nWat + 2 * $nCLA)
}

proc ionmolar2 { mol {resname CLA} } {
    # Finds molar concentration of ions
    set selWat [atomselect $mol "name OH2"]
    set nWat [$selWat num]
    $selWat delete
    
    set selCLA [atomselect $mol "resname $resname"]
    set nCLA [$selCLA num]
    $selCLA delete
    
    puts "N(water) = $nWat"
    puts "N($resname) = $nCLA"
    
    expr 55.6 * $nCLA / ($nWat + $nCLA)
}

proc distributeMols { args } {
    set molids	$args
    set num	[llength $args]
    
    # Parameters
    set xmax 1.85	;# These could be more useful numbers ...
    set ymax 1.5
    
    # Translate and fix
    set i 0
    foreach molid $molids {
	set dx [expr (($i+0.5)/$num) * 2 * $xmax - $xmax]
	puts "dx = $dx"
	translate to $dx 0 0
	mol fix $molid
	incr i
    }
    
    # Finally, unfix all
    foreach molid $molids {
	mol free $molid
    }
}

proc clip { clipid repid molid func arg } {
    mol clipplane $func $clipid $repid $molid $arg
}


proc canon { filename } {
    return [exec readlink -f $filename]
}

### STACK PROCS ###
#
proc push { stack value } {
    upvar $stack list
    lappend list $value
}

proc pop { stack } {
    upvar $stack list
    set value [lindex $list end]
    set list [lrange $list 0 [expr [llength $list]-2]]
    return $value
}

proc shift { stack } {
    upvar $stack list
    set value [lindex $list 0]
    set list [lrange $list 1 end]
    return $value
}

proc unshift { stack value } {
    upvar $stack list
    set list [concat $value $list]
}

proc stackrotate { stack {value 1} } {
    upvar $stack list
    while { $value > 0 } {
	set el [shift list]
	push list $el
	incr value -1
    }
    while { $value < 0 } {
	set el [pop list]
	unshift list $el
	incr value 1
    }
    return $list
}

proc ul { str } {
    # This proc brackets the provided string with the terminal escape sequence for underlining
    return "\033\[4m$str\033\[0m"
}

proc silent { } { return }

proc vmdargs { script args } {
    upvar argv argv_local
    
    if { [llength $argv_local] != [llength $args] } {
	puts -nonewline "Usage: vmd -dispdev text -e $script -args"
	foreach arg $args {
	    puts -nonewline " [ul $arg]"
	}
	puts ""
	exit -1
    }
    
    foreach arg $args {
	upvar $arg arg_local
	set arg_local [shift argv_local]
    }
}

proc vmdargslist { script args } {
    upvar argv argv_local
    
    # DEBUG
    puts "argv: $argv_local"
    
    set nargs [llength $args]
    set nlistargs 0
    set i 0
    foreach arg $args {
	# Loop through and check for a list arg ... more than one is a no-no
	if { [string match "@*" $arg] } {
	    incr nlistargs
	    set listargind $i
	}
	incr i
    }
    if { $nlistargs > 1 } {
	puts "Cannot use more than one list variable in vmdargs call!"
	exit -1
    }
    
    if { ($nlistargs == 0 && [llength $argv_local] != [llength $args]) || [llength $argv_local] < [llength $args] } {
	puts -nonewline "Usage: vmd -dispdev text -e $script -args"
	foreach arg $args {
	    if { [string match "@*" $arg] } {
		set arg2 [string range $arg 1 end]
		puts -nonewline " [ul $arg2] \[ [ul $arg2] ... \]"
	    } else {
		puts -nonewline " [ul $arg]"
	    }
	}
	puts ""
	exit -1
    }
    
    foreach arg $args {
	if { [string match "@*" $arg] } {
	    # remove leading '@' symbol
	    set arg2 [string range $arg 1 end]
	    upvar $arg2 arg_local
	    while { [llength $argv_local] >= $nargs - $listargind } {
		lappend arg_local [shift argv_local]
	    }
	} else {
	    upvar $arg arg_local
	    set arg_local [shift argv_local]
	}
    }
}

proc vmdusage { script args } {
    puts -nonewline "Usage: vmd -dispdev text -e $script -args"
    foreach arg $args {
	puts -nonewline " [ul $arg]"
    }
    puts ""
}

proc s2hms { s } {
    set h [expr $s/3600]
    set m [format %02d [expr $s % 3600 / 60]]
    set s [format %02d [expr $s % 60]]
    return "$h:$m:$s"
}

set progressbar_t0 0
set progressbar_t1 0
proc progressbar { i n {length 50} } {
    global progressbar_t0 progressbar_t1
    
    if { $length == 0 } { set length -1 }
    
    if { $length > 0 } {
	set ii [expr $i*$length/$n]	;# rescale variables
	set iip1 [expr ($i+1)*$length/$n]
	set nn $length
    } else {
	set ii [expr -$i/$length]
	set iip1 [expr -($i+1)/$length]
    }
    
    if { $i == 0 } {
	set progressbar_t0	[clock seconds]
	set t			0
	set progressbar_t1	0
    } else {
	set t [expr [clock seconds] - $progressbar_t0]
    }
    
    if { $iip1 == $ii && $i != 0 && $i != $n-1 && ($t == $progressbar_t1 || $length < 0) } {
	# Don't draw if this would draw the same thing as the last iteration
	return
    }
    
    if { $length > 0 } {
	puts -nonewline "\r|"
	for { set j 0 } { $j < $nn } { incr j } {
	    if { $j < $iip1 } {
		puts -nonewline "="
	    } else {
		puts -nonewline " "
	    }
	}
    }
    
    if { $i == 0 } {
	set time "99:99:99"
    } elseif { $i == [expr $n-1] } {
	set time [s2hms $t]
    } else {
	# Calculate estimated remaining time
	set rt [expr $t * ($n-$i) / $i]
	set time [s2hms $rt]
    }
    
    if { $length > 0 } {
	puts -nonewline "|  [expr 100*($i+1)/$n]%  [expr $i+1] / $n  $time "
	flush stdout
	if { $i == [expr $n-1] } {
	    puts ""
	}
    } else {
	puts "[expr $i+1] / $n  $time"
    }
    
    set progressbar_t1 $t
}

proc vecadd_list { veclist1 veclist2 } {
    set l1 [llength $veclist1]
    set l2 [llength $veclist2]
    
    if { $l1 == 0 } {
	return $veclist2
    } elseif { $l2 == 0 } {
	return $veclist1
    } elseif { $l1 != $l2 } {
	puts "Lists must have same length! ($l1 != $l2)"
	return -1
    }
    
    set veclist3 {}
    foreach vec1 $veclist1 vec2 $veclist2 {
	lappend veclist3 [vecadd $vec1 $vec2]
    }
    
    return $veclist3
}

proc vecsub_list { veclist1 veclist2 } {
    set l1 [llength $veclist1]
    set l2 [llength $veclist2]
    
    set veclist3 {}
    
    if { $l1 == 0 } {
	return $veclist2
    } elseif { $l2 == 0 } {
	return $veclist1
    } elseif { $l2 == 3 } {
	set vec2 $veclist2
	foreach vec1 $veclist1 {
	    lappend veclist3 [vecsub $vec1 $vec2]
	}
	return $veclist3
    } elseif { $l1 != $l2 } {
	puts "Lists must have same length! ([llength $veclist1] != [llength $veclist2])"
	return -1
    } else {
	foreach vec1 $veclist1 vec2 $veclist2 {
	    lappend veclist3 [vecsub $vec1 $vec2]
	}
	return $veclist3
    }
}

proc list_sqrt { list1 } {
    set results {}
    foreach el $list1 {
	lappend results [expr sqrt($el)]
    }
    return $results
}

proc list_add { list1 delta } {
    set results {}
    foreach el $list1 {
	lappend results [expr {$el + $delta}]
    }
    return $results
}

proc inner_product { list1 list2 } {
    if { [llength $list1] != [llength $list2] } {
	return "nan"
    }
    set sum 0.
    foreach el1 $list1 el2 $list2 {
	set sum [expr {$sum + $el1 * $el2}]
    }
    return $sum
}

proc lscale { list1 scale } {
    set results {}
    foreach el $list1 {
	lappend results [expr $el * $scale]
    }
    return $results
}

proc lapply { proc list { arg1 "" } } {
    set results {}
    switch [llength $arg1] \
	0 { foreach el $list { lappend results [$proc $el] } } \
	1 { foreach el $list { lappend results [$proc $el $arg1] } } \
	"[llength $list]" { foreach el1 $list el2 $arg1 { lappend results [$proc $el1 $el2] } } \
	default { return -1 }
    return $results
}

proc vecscale_list { veclist scale } {
    return [lapply vecscale $veclist [expr $scale]]
}

proc veclength_list { veclist } {
    return [lapply veclength $veclist]
}

proc veclength2_list { veclist } {
    return [lapply veclength2 $veclist]
}

proc vecsum_list { veclist } {
    set sum "0. 0. 0."
    foreach vec $veclist {
	set sum [vecsub $sum $vec]
    }
    return $sum
}

proc vecmean_list { veclist } {
    set sum [vecsum_list $veclist]
    set mean [vecscale $sum [expr {1./[llength $veclist]}]]
    return $mean
}

proc coordtrans_list { mat coordlist1 } {
    set coordlist2 {}
    foreach coord $coordlist1 {
	lappend coordlist2 [coordtrans $mat $coord]
    }
    return $coordlist2
}

proc lrotate { list1 } {
    # e.g. rotate_list "a b c" -> "c a b"
    upvar $list1 list
    set a [lrange $list 0 [expr [llength $list]-2]]
    set b [list [lindex $list end]]
    set list [concat $b $a]
}

proc lsum list {
    # from wiki.tcl.tk/951
    # UPDATE this method totally sucks
    #expr [join $list +] + 0
    set sum 0.0
    foreach el $list {
	set sum [expr $sum + $el]
    }
    return $sum
}

# proc lmean list {
#     if { [llength $list] == 0 } {
# 	return 0
#     }
#     expr ([join $list +])/[llength $list].0
# }
proc lmean list {
    if { [llength $list] == 0 } {
	return 0
    }
    set sum 0.0
    foreach el $list {
	set sum [expr {$sum + $el}]
    }
    return [expr {$sum/[llength $list]}]
}

proc lmean_w { list weights } {
    if { [llength $list] == 0 } {
	return 0
    }
    if { [llength $list] != [llength $weights] } {
	return 0
    }
    set sum 0.0
    foreach el $list w $weights {
	set sum [expr $sum + $el*$w]
    }
    return [expr $sum/[lsum $weights]]
}

proc lstddev list {
    if { [llength $list] < 2 } {
	return 0
    }
    
    set mean [lmean $list]
    set var 0.0
    foreach el $list {
	set var [expr {$var + (($el-$mean)*($el-$mean))}]
    }
    
    return [expr {sqrt($var/[expr [llength $list]-1])}]
}

proc lstderr list {
    if { [llength $list] < 2 } {
	return 0
    }
    
    return [expr {[lstddev $list]/sqrt([expr [llength $list]-1])}]
}

proc lstddev_w { list weights } {
    if { [llength $list] < 2 } {
	return 0
    }
    if { [llength $list] != [llength $weights] } {
	return 0
    }
    
    set mean [lmean_w $list $weights]
    set var 0.0
    foreach el $list w $weights {
	set var [expr $var + ($w*($el-$mean)*($el-$mean))]
    }
    
    return [expr sqrt($var / ([llength $list]-1) / [lmean $weights])]
}

proc lmin { list } {
    lindex [lsort -real $list] 0
}

proc lmax { list } {
    lindex [lsort -real $list] end
}


# proc lintersect args {
#     # From http://wiki.tcl.tk/43
#     set res {}
#     foreach element [lindex $args 0] {
# 	set found 1
# 	foreach list [lrange $args 1 end] {
# 	    if {[lsearch -exact $list $element] < 0} {
# 		set found 0; break
# 	    }
# 	}
# 	if {$found} {lappend res $element}
#     }
#     set res
# } ;# RS
proc lintersect args {
    set res {}
    foreach element [lindex $args 0] {
	if { $element in [lindex $args 1] } {
	    lappend res $element
	}
    }
    return $res
}


proc trajdev { ref sel { frame0 0 } } {
    # Sets user field to deviation of sel2 from sel1, after first aligning the structures
    set molid [$sel molid]
    set nframes [molinfo $molid get numframes]
    set selall [atomselect $molid all]
    for { set frame $frame0 } { $frame < $nframes } { incr frame } {
	progressbar $frame $nframes -10
	$sel frame $frame
	
	# Align sel with ref
	set mat [measure fit $sel $ref]
	$selall frame $frame
	$selall move $mat
	
	# Find (root square) deviation
	set xref [$ref get {x y z}]
	set xsel [$sel get {x y z}]
	set dx [lapply vecsub $xref $xsel]
	set dev [lapply veclength $dx]
	
	# Set user field
	$sel set user $dev
    }
}

proc trajdev2 { ref sel { frame0 0 } } {
    # Sets user field to SQUARE deviation of sel2 from sel1, after first aligning the structures
    set molid [$sel molid]
    set nframes [molinfo $molid get numframes]
    set selall [atomselect $molid all]
    for { set frame $frame0 } { $frame < $nframes } { incr frame } {
	progressbar $frame $nframes -10
	$sel frame $frame
	
	# Align sel with ref
	set mat [measure fit $sel $ref]
	$selall frame $frame
	$selall move $mat
	
	# Find (root square) deviation
	set xref [$ref get {x y z}]
	set xsel [$sel get {x y z}]
	set dx [lapply vecsub $xref $xsel]
	set dev2 [lapply veclength2 $dx]
	
	# Set user field
	$sel set user $dev2
    }
}

proc trajdev_noalign { ref sel { frame0 0 } } {
    # Sets user field to deviation of sel2 from sel1, without aligning the structures
    set molid [$sel molid]
    set nframes [molinfo $molid get numframes]
    set selall [atomselect $molid all]
    for { set frame $frame0 } { $frame < $nframes } { incr frame } {
	progressbar $frame $nframes -10
	$sel frame $frame
	
	# Find (root square) deviation
	set xref [$ref get {x y z}]
	set xsel [$sel get {x y z}]
	set dx [lapply vecsub $xref $xsel]
	set dev [lapply veclength $dx]
	
	# Set user field
	$sel set user $dev
    }
}


### NAMD-RELATED PROCEDURES ###
#
proc dcdFrames { dcdfile } {
    return [exec catdcd -num $dcdfile | awk "/^Total frames:/ \{ print \$3 \}"]
}

proc get_first_ts { xscfile } {
    set fd [open $xscfile r]
    gets $fd
    gets $fd
    gets $fd line
    set ts [lindex $line 0]
    close $fd
    return $ts
}

proc readExtendedSystem { xsfile } {
    if { [string match *.xst $xsfile] } {
	set xst 1
    } elseif { [string match *.xsc $xsfile] } {
	set xst 0
    } else {
	return -1
    }
    
    set ch [open $xsfile r]
    
    set xsdata {}
    while { [gets $ch line] != -1 } {
	if { [regexp "^#" $line] } { continue }
	if { ! $xst } {
	    return [split $line]
	} else {
	    lappend xsdata [split $line]
	}
    }
    return $xsdata
}

# like readExtendedSystem, but puts data in array
proc readExtendedSystem_array { xsfile arrayvar } {
    upvar $arrayvar xsdata
    
    if { [string match *.xst $xsfile] } {
	set xst 1
    } elseif { [string match *.xsc $xsfile] } {
	set xst 0
    } else {
	return -1
    }
    
    set ch [open $xsfile r]
    
    #set xsdata {}
    while { [gets $ch line] != -1 } {
	if { [regexp "^#" $line] } { continue }
	if { ! $xst } {
	    return [split $line]
	} else {
	    set cols [split $line]
	    set xsdata([lindex $cols 0]) [lrange $cols 1 end]
	}
    }
    #return $xsdata
}


### FROM JEFF ###
#
proc writeExtendedSystem {xscList outFile} {
    set out [open $outFile w]
    puts $out "# NAMD extended system configuration restart file"
    puts $out "#\$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w"
    puts $out $xscList
    close $out
    return
}

# Get the position of an atom.
proc getPos {segName resId name {mole top}} {
    set sel [atomselect $mole "segname $segName and resid $resId and name $name"]
    set n [$sel num]
    
    if {$n < 1} {
	puts "Warning! Atom ${segName}:${resId}:${name} does not exist."
	return [list 0.0 0.0 0.0]
    } elseif {$n > 1} {
	puts "Warning! Atom ${segName}:${resId}:${name} in not unique."
    }

    set r [lindex [$sel get {x y z}] 0]
    $sel delete
    return $r
}

# Get the standard position of a DNA basepair.
proc getBasepairPos {segA resA segB resB {mole top}} {
    set car1A [getPos $segA $resA "C1'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    return [vecScale 0.5 [vecAdd $car1A $car1B]] 
}

# Define a basis for a base.
proc getBaseNormal {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
    # Get the residue name.
    set sel [atomselect top $selText]
    set resName [lindex [$sel get resname] 0]
    $sel delete

    # Get the hexagon ring basis.
    if {[string equal ADE $resName] || [string equal GUA $resName]} {
	set selX0 [atomselect $mole "($selText) and name C4"]
	set selX1 [atomselect $mole "($selText) and name N1"]
	set selY0 [atomselect $mole "($selText) and name N3"]
	set selY1 [atomselect $mole "($selText) and name C5"]
    } else {
	set selX0 [atomselect $mole "($selText) and name N1"]
	set selX1 [atomselect $mole "($selText) and name C4"]
	set selY0 [atomselect $mole "($selText) and name C2"]
	set selY1 [atomselect $mole "($selText) and name C6"]
    }

    set rX0 [lindex [$selX0 get {x y z}] 0]
    set rX1 [lindex [$selX1 get {x y z}] 0]
    set rY0 [lindex [$selY0 get {x y z}] 0]
    set rY1 [lindex [$selY1 get {x y z}] 0]

    $selX0 delete
    $selX1 delete
    $selY0 delete
    $selY1 delete

    set ex [vecsub $rX1 $rX0]
    set ex [vecscale [expr 1.0/[veclength $ex]] $ex]
    set ey [vecsub $rY1 $rY0]
    set ey [vecsub $ey [vecscale [vecdot $ey $ex] $ex]]
    set ey [vecscale [expr 1.0/[veclength $ey]] $ey]
    set ez [veccross $ex $ey]

    return $ez
}


# Define a basis for a base.
proc getBasepairDirection {segA resA segB resB {mole top}} {
    set car1A [getPos $segA $resA "C1'" $mole]
    set car1B [getPos $segB $resB "C1'" $mole]

    return [vecUnit [vecsub $car1B $car1A]]
}

# Define a basis for a basepair.
proc getBasepairBasisBest {segA resA segB resB {mole top}} {
    set zA [getBaseNormal $segA $resA $mole]
    set zB [getBaseNormal $segB $resB $mole]
    set ez [vecUnit [vecsub $zA $zB]]

    set ex [getBasepairDirection $segA $resA $segB $resB $mole]
    set ex [vecUnit [vecsub $ex [vecscale [vecdot $ex $ez] $ez]]]
    
    set ey [veccross $ez $ex]
    return [matTranspose [list $ex $ey $ez]]
}


# Get the standard position of DNA nucleotide.
proc getNucleotidePos {seg res {mole top}} {
    set rC5P [getPos $seg $res "C5'" $mole]
    set rC3P [getPos $seg $res "C3'" $mole]

    return [vecScale 0.5 [vecAdd $rC3P $rC5P]] 
}

# Define a basis for the nucleotide.
proc getNucleotideBasis {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
    set selO5P [atomselect $mole "($selText) and name O5'"]
    set selO3P [atomselect $mole "($selText) and name O3'"]
    set selC1P [atomselect $mole "($selText) and name C1'"]

    set resName [lindex [$selO5P get resname] 0]
    if {[string equal ADE $resName] || [string equal GUA $resName]} {
	set selN [atomselect $mole "($selText) and name N9"]
    } else {
	set selN [atomselect $mole "($selText) and name N1"]
    }

    set rO5P [lindex [$selO5P get {x y z}] 0]
    set rO3P [lindex [$selO3P get {x y z}] 0]
    set rC1P [lindex [$selC1P get {x y z}] 0]
    set rN [lindex [$selN get {x y z}] 0]
    
    $selO5P delete
    $selO3P delete
    $selC1P delete
    $selN delete

    set ex [vecsub $rO3P $rO5P]
    set ex [vecscale [expr 1.0/[veclength $ex]] $ex]
    set ey [vecsub $rC1P $rO5P]
    set ey [vecsub $ey [vecscale [vecdot $ey $ex] $ex]]
    set ey [vecscale [expr 1.0/[veclength $ey]] $ey]
    set ez [veccross $ex $ey]

    return [matTranspose [list $ex $ey $ez]]
}

# Get the standard position of DNA base.
proc getBasePos {seg res {mole top}} {
    return [getPos $seg $res "C1'" $mole]
}

# Get the center of mass of DNA base
proc getBaseCOM {seg res {mole top}} {
    set sel [atomselect $mole "segname $seg and resid $res and nucleic"]
    set com [measure center $sel weight mass]
    $sel delete
    return $com
}

# Define a basis for the base.
proc getBaseBasis {segName resId {mole top}} {
    set selText "segname $segName and resid $resId"
    set sel [atomselect $mole $selText]
    set resName [lindex [$sel get resname] 0]
    $sel delete

    # Get the hexagon ring basis.
    if {[string equal ADE $resName] || [string equal GUA $resName]} {
	set selX0 [atomselect $mole "($selText) and name C4"]
	set selX1 [atomselect $mole "($selText) and name N1"]
	set selY0 [atomselect $mole "($selText) and name C4"]
	set selY1 [atomselect $mole "($selText) and name C6"]
    } else {
	set selX0 [atomselect $mole "($selText) and name C6"]
	set selX1 [atomselect $mole "($selText) and name N3"]
	set selY0 [atomselect $mole "($selText) and name C6"]
	set selY1 [atomselect $mole "($selText) and name C4"]
    }

    set rX0 [lindex [$selX0 get {x y z}] 0]
    set rX1 [lindex [$selX1 get {x y z}] 0]
    set rY0 [lindex [$selY0 get {x y z}] 0]
    set rY1 [lindex [$selY1 get {x y z}] 0]

    $selX0 delete
    $selX1 delete
    $selY0 delete
    $selY1 delete
      
    set ex [vecsub $rX1 $rX0]
    set ex [vecscale [expr 1.0/[veclength $ex]] $ex]
    set ey [vecsub $rY1 $rY0]
    set ey [vecsub $ey [vecscale [vecdot $ey $ex] $ex]]
    set ey [vecscale [expr 1.0/[veclength $ey]] $ey]
    set ez [veccross $ex $ey]

    return [matTranspose [list $ex $ey $ez]]
}


### PROCS WHICH NAMD HAS DEFINED (FOR TCLFORCES/TCLBC) ###
#
proc atomid { segname resid atomname } {
    set sel [atomselect top "segname $segname and resid $resid and name $atomname"]
    set index [$sel get index]
    $sel delete
    return $index
}

proc addatom { atomid } {
    return
}

proc print { text } {
    puts $text
}

