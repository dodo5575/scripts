## Plot data (as fn of time) in a vmd script scripts

proc check_molid { ID } {
    if { [string equal $ID top] } { return 1 }
    foreach MID [molinfo list] {
	if { $ID == $MID } { return 1 }
    }
    error "MolID $ID does not exist" 
}
proc check_molID { ID } { return [check_molid $ID] }	 

proc shift { listName } {
    ## like bash cmd "shift"
    upvar $listName list
    set list [lrange $list 1 end]
}

namespace eval tplot {
    variable subspace_count 0
    variable subspaces "" ;# list of active (sub)namespaces
    
    proc new { args } {
	## Handle Args
	set args [join $args]

	set dataFile [lindex $args 0]
	shift args
	
	set opts "{mol top} {column end} {color Red}"
	foreach {var val} [join $opts] { set $var $val }	    
	foreach {var val} $args { set $var $val }
		
	if { [string equal $mol top] } { set mol [molinfo top] }
	check_molID $mol
	set molID $mol
	
	## create a new namespace for the plot
	variable subspaces
	variable subspace_count
	variable subspace tplot$subspace_count
	lappend subspaces $subspace
	incr subspace_count
	
	## create a new namespace that contains data for the plot
	namespace eval $subspace {
	    puts [namespace current]
	    variable data ""
	    variable max ""
	    variable min ""
	    variable drawID ""
	    variable highlight ""
	    
	    ## define things that user can get/set
	    variable width 1.25;
	    variable height 1.;
	    variable x .3;
	    variable y -1.;
	    variable z .97;
	    variable r .015;
	    
	    variable timePerFrame .05
	    
 	    variable xtitle "time (ns)"
	    variable ytitle ""
	    
	    ## create graphics molecule
	    variable drawID [mol load graphics graphics]
	}
	set ${subspace}::molID $molID
	mol top $molID
	
	## reset view/fix draw molecule
	variable ${subspace}::drawID
	resetView $subspace
	mol inactive $drawID
	mol fix $drawID


	## create proc in global namespace that links to user interface
	proc ::$subspace { args } "::tplot::subspace $subspace \[join \$args\]"
	
	## gather data
	load_data $subspace $dataFile $column
	
	## plot static data
	draw $subspace
	#doTrace $subspace
	
	## highlight data point
	## add a dummy traceproc
	proc ${subspace}::doTrace { array ID op } "::tplot::doTrace $subspace \$array \$ID \$op"
	#puts "trace variable ::vmd_frame($molID) w ${subspace}::doTrace"
	trace variable ::vmd_frame($molID) w ::tplot::${subspace}::doTrace 

	## draw
	subspace $subspace draw

	## return the namespace name 
	return $subspace

    } ;## end new

    ## USER INTERFACE!!!
    proc subspace { subspace args } {
	set args [join $args]
	#puts "args are $args"
	## define possible commands
	set rm "rm|remove|delete"
	set draw "draw|refresh|redraw"
	
	## get command
	set cmd [lindex $args 0] 
	switch $cmd "$rm { set cmd rm }"
	set args [lrange $args 1 end]
	## execute command
	switch $cmd {
	    rm { ::tplot::rm $subspace $args }
	    moveto { ::tplot::moveto $subspace $args }
	    draw { ::tplot::draw $subspace $args }
	    
	    default {
		puts "tplot: Unknown command $cmd.  Valid commands are $rm & $draw."
	    }
	}
    }


    
    ## "private" procs (can't actually make things private)
    proc load_data { subspace files {column end} } {
	puts "column is $column"
	## read data in $file to data
	variable ${subspace}::data ""
	foreach file $files {
	    if { ![file exists $file] } {
		puts "WARNING: file $file does not exist... skipping"
		continue
	    }
	    set ch [open $file r]
	    while { [gets $ch line] > 0 } {
		if { [regexp "^#" $line] } {
		    continue
		}
		set line [split $line]
		lappend data [lindex $line $column]
	    }
	    close $ch
	}
	
	## average data to fit to numframes
	set dl [llength $data]
	set fl [molinfo top get numframes]
	set data2 ""
	if { $dl > $fl } {
	    puts "WARNING: data is $dl long and does not match $fl frames... averaging"
	    for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
		lappend data2 [vecmean [lrange $data [expr round($i*$dl/$fl)] [expr round(($i+1)*$dl/$fl)]]]
	    }
	    set data $data2
	}
	## define lmax/lmin
	variable ${subspace}::max [::lmax $data]
	variable ${subspace}::min [::lmin $data]
	
	#puts "minmax were $min $max"
	foreach {min max div} [getRange $min $max] { break }
	#puts "minmax are $min $max"

	## scale data
	set data2 ""
	foreach d $data {
	    lappend data2 [expr ($d-$min)/($max-$min)]
	}
	set data $data2
    }
    proc getRange {min max} {
	set range [expr $max-$min]
	set pow [expr floor(log10($range))-1]
	set range [expr $range / (10**$pow)]

	set divmin [expr $range/5]
	set div [expr min( (1 > $divmin ? 1 : 100), (2.5 > $divmin ? 2.5 : 100),  (5 > $divmin ? 5 : 100), (10 > $divmin ? 10 : 100), (25 > $divmin ? 25 : 100) )]
	
	set range [expr $range * 10**$pow]
	set div  [expr $div * 10**$pow]
	
	set min [expr floor($min/$div)*$div]
	set max [expr  ceil($max/$div)*$div]
	return "$min $max $div"
    }
    proc resetView { subspace } {
	variable ${subspace}::drawID
	foreach matrix {center_matrix scale_matrix rotate_matrix global_matrix} {
	    molinfo $drawID set $matrix "{[transidentity]}"
	}
    }
	
    proc draw { subspace args } {
	foreach var "drawID min max data width height x y z xtitle ytitle timePerFrame" { variable ${subspace}::$var }
	graphics $drawID delete all
	## draw plot background
	graphics $drawID color white
	graphics $drawID triangle "$x $y $z" "[expr {$x+$width}] $y $z" "[expr {$x+$width}] [expr {$y+$height}] $z"
	graphics $drawID triangle "$x $y $z" "$x [expr {$y+$height}] $z" "[expr {$x+$width}] [expr {$y+$height}] $z"

	## draw axes
	graphics $drawID color black
	graphics $drawID line "$x $y $z" "[expr {$x+$width}] $y $z"
	graphics $drawID line "[expr {$x+$width}] $y $z" "[expr {$x+$width}] [expr {$y+$height}] $z"
	graphics $drawID line "[expr {$x+$width}] [expr {$y+$height}] $z" "$x [expr {$y+$height}] $z"
	graphics $drawID line "$x [expr {$y+$height}] $z" "$x $y $z" 
	
	## draw tics
	## draw data
	set z [expr $z+.01]
	graphics $drawID color red
	set i 0
	set imax [llength $data]
	set d1 [lindex $data 0]
	for {set i 1} {$i < [llength $data]} {incr i} {
	    set d2 [lindex $data $i]
	    graphics $drawID line "[expr {$x + ($i-1)*$width/$imax}] [expr {$y + $d1*$height}] $z"  "[expr {$x + $i*$width/$imax}] [expr {$y + $d2*$height}] $z" 
	    set d1 $d2
	}
	## draw labels
	graphics $drawID color black
	set c_width .056 ;# works well with whatever my default is
	set c_width .04 ;# works well with display resize 800 800
	
	set s [format "%g" 0]
	set l [string length $s]
	graphics $drawID text "[expr $x+$width*(0-.02*$l)] [expr $y-$c_width*$height] $z" $s
	set s [format "%g" [expr [llength $data]*$timePerFrame]]
	set l [string length $s]
	graphics $drawID text "[expr $x+$width*(1-.02*$l)] [expr $y-$c_width*$height] $z" $s
	
	set s [format "%g" $min]
	set l [string length $s]
	graphics $drawID text "[expr $x-$width*($c_width*($l+.5))] $y $z" $s
	set s [format "%g" $max]
	set l [string length $s]
	graphics $drawID text "[expr $x-$width*($c_width*($l+.5))] [expr $y+$height] $z" $s
	
	## draw title
	set s $xtitle
	set l [string length $s]
	graphics $drawID text "[expr $x+$width*(.5-.5*$c_width*$l)] [expr $y- 2*$c_width*$height] $z" $s
	set s $ytitle
	set l [string length $s]
	graphics $drawID text "[expr $x+0*$width] [expr $y+ (1+1.5*$c_width)*$height] $z" $s
	set z [expr $z-.01]
	
	doTrace $subspace
    }
    proc doTrace { subspace args } {
	## this could be cleaner...
	foreach var "molID drawID highlight min max data width height x y z r" { variable ${subspace}::$var }
	
	set z [expr $z + .02]
	set i $::vmd_frame($molID)
	set imax [llength $data]
	set d1 [lindex $data $i]
	## update trajectory plot each frame
	graphics $drawID color blue
	set replace_code ""
	if { $highlight != "" } {
	    foreach h $highlight {
		lappend replace_code "graphics $drawID replace $h"
	    }
	}
	eval [lindex $replace_code 0]
	lappend highlight [graphics $drawID sphere "[expr {$x + $i*$width/$imax}] [expr {$y + $d1*$height}] $z" radius $r resolution 16]
	eval [lindex $replace_code 1] 
	lappend highlight [graphics $drawID line "[expr {$x + $i*$width/$imax}] $y $z" "[expr {$x + $i*$width/$imax}] [expr $y+$height] $z"]
	eval [lindex $replace_code 2] 
	lappend highlight [graphics $drawID line "$x [expr {$y+$d1*$height}] $z" "[expr {$x + $width}] [expr {$y+$d1*$height}] $z"]
	# 	    eval [lindex $replace_code 3] 
	#	    lappend highlight [graphics $drawID text "[expr {$x+$width*.9}] [expr {$y+1.1*$height}] $z" [format "%.4g" [expr $d1*$max - $min]]]
	set z [expr $z - .02]
    }
    proc rm { subspace args } {
	foreach var "drawID trace" { variable ${subspace}::$var }
		    
	## remove plot molecule
	mol delete $drawID
	
	## remove trace_plot
	
	## remove namespace
	namespace delete $subspace
    }
    proc moveto { subspace args } {
	foreach {x y} [join $args] { break }
	set ${subspace}::x [expr {double($x)}]
	set ${subspace}::y [expr {double($y)}]
	graphics [subst $${subspace}::drawID] delete all
	draw $subspace
	doTrace $subspace
    }
}

## Direct user interface
proc tplot { args } {
    ## first arg is a command
    set new "new|plot"
    
    set cmd [lindex $args 0] 
    puts "cmd is $cmd"
    puts "tplot.args is $args"
    switch $cmd "$new { set cmd new }"
    #puts "tplot::new [lindex $args 1] [lrange $args 2 end]"
    switch $cmd {
	new { tplot::new [lrange $args 1 end] }
	default {
	    puts "tplot: Unknown command [lindex $args 0].  Valid commands are $new"
	}
    }
}