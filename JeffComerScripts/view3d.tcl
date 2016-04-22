# This script implements three-dimesional view.
# Author: Jeff Comer <jcomer2@illinois.edu>
package provide view3d 1.0

proc view3d {} {
	::View3d::view3dGUI	
}

namespace eval ::View3d:: {
	namespace export view3dGUI
	
	variable w
	variable basis
	variable psi
	variable mouse_x
	variable mouse_y
	
	variable points
	variable colors
	variable size
	variable zoom
	variable z0
	variable dz
	
	variable axes
	variable acolors
	variable asize
	variable azoom
	variable az0
	variable adz
}

proc ::View3d::initialize {} {
	variable w
	variable basis
	variable psi
	variable mouse_x
	variable mouse_y
	
	variable points
	variable colors
	variable size
	variable zoom
	variable z0
	variable dz
	
	variable axes
	variable acolors
	variable asize
	variable azoom
	variable az0
	variable adz
		
	set basis [transaxis x 1.0 pi]
	set psi [expr 4.0*atan(1.0)/100.0]
	set mouse_x 0.0
	set mouse_y 0.0
	
	set points [list {1 1 1} {1 1 -1} {1 -1 1} \
		{1 -1 -1} {-1 1 1} {-1 1 -1} {-1 -1 1} {-1 -1 -1}]
	lappend points {1 0 0} {-1 0 0} {0 1 0} {0 -1 0} {0 0 1} {0 0 -1}
	lappend points {0 0 0}
	for {set i 0} {$i < 20} {incr i} {
		lappend points \
		[list [expr rand()-0.5] [expr rand()-0.5] [expr rand()-0.5]]
	}
	
	set colors {}
	set size 320
	set zoom [expr 0.2*$size]
	set z0 [expr -[range $points]]
	set dz [expr -2.0*$z0]
	
	foreach r $points {lappend colors {0 0 0}}	
	
	set axes [list {1 0 0} {2 0 0} {3 0 0} {4 0 0} {5 0 0} \
		{0 1 0} {0 2 0} {0 3 0} {0 4 0} {0 5 0} \
		{0 0 1} {0 0 2} {0 0 3} {0 0 4} {0 0 5} {0 0 0}]
	set acolors [list {1 0 0} {1 0 0} {1 0 0} {1 0 0} {1 0 0} \
		{0 1 0} {0 1 0} {0 1 0} {0 1 0} {0 1 0} \
		{0 0 1} {0 0 1} {0 0 1} {0 0 1} {0 0 1} {0 0.5 0.5}]
	set asize 20
	set azoom [expr 0.2*$asize]
	set az0 [expr -[range $axes]]
	set adz [expr -2.0*$az0]
}

proc ::View3d::view3dGUI {} {
	variable w
	variable size

	initialize
	if { [winfo exists .view3dgui] } {
		wm deiconify .view3dgui
		return
	}
	
	set w [toplevel ".view3dgui"]
	wm title $w "View3d"
	
	frame $w.mbar -borderwidth 1 -relief raised
	pack $w.mbar -fill x
	menubutton $w.mbar.file -text "File" -menu $w.mbar.file.m
	pack $w.mbar.file -side left
	menu $w.mbar.file.m
	$w.mbar.file.m add command -label "Open"

	menubutton $w.mbar.edit -text "Edit" -menu $w.mbar.edit.m
	pack $w.mbar.edit -side left
	menu $w.mbar.edit.m
	$w.mbar.edit.m add command -label "Clear" -command ::View3d::clear

	frame $w.style -borderwidth 1 -relief sunken
	pack $w.style -fill x
	
	label $w.style.readout -text "x: 0.00 y: 0.00"
	pack $w.style.readout -side right

	canvas $w.viewport -background white -width $size -height $size
	bind $w.viewport <ButtonPress-1> {::View3d::track %x %y}
	bind $w.viewport <B1-Motion> {::View3d::drag %x %y}
	pack $w.viewport
	redraw
	
	return $w
}

# Take colors (r, g, b) on the range [0,1) and return the
# color in the #rrggbb hexadecimal format.
proc ::View3d::getColor {r g b} {
	set r [expr int(255.0*$r)]
	set g [expr int(255.0*$g)]
	set b [expr int(255.0*$b)]
	
	if {$r < 0} {set r 0}
	if {$g < 0} {set g 0}
	if {$b < 0} {set b 0}
	if {$r > 255} {set r 255}
	if {$g > 255} {set g 255}
	if {$b > 255} {set b 255}
	
	set sr ""
	set sg ""
	set sb ""
	if {$r < 16} {set sr 0}
	if {$g < 16} {set sg 0}
	if {$b < 16} {set sb 0}
	return [format "#%s%x%s%x%s%x" $sr $r $sg $g $sb $b]
}

proc ::View3d::range {pt} {
	set dmax 0.0
	foreach r $pt {
		set d [veclength $r]
		if {$d > $dmax} {set dmax $d}
	}
	return $dmax
}

proc ::View3d::render {points colors cx cy zoom z0 dz} {
	variable w
	variable basis	
	
	# Project the points orthographically.
	set pts {}
	foreach r $points c $colors {
		lappend pts [concat [vectrans $basis $r] $c]
	}
	
	# Sort by distance and determine the colors.
	set pts [lsort -real -decreasing -index 2 $pts]
				
	# Draw the points.
	foreach p $pts {
		foreach {x y z r g b} $p {break}
		set sx [expr $x*$zoom + $cx]
		set sy [expr $y*$zoom + $cy]
		set in [expr 0.9*($z-$z0)/$dz]
		set c1 [getColor [expr $r+$in] [expr $g+$in] [expr $b+$in]]
		$w.viewport create rectangle [expr $sx-2] [expr $sy-2] \
			[expr $sx+2] [expr $sy+2] -fill $c1 -outline $c1
	}
}

proc ::View3d::track {mx my} {
	variable w
	variable mouse_x
	variable mouse_y
	
	#Set the new coordinates.
	set mouse_x $mx
	set mouse_y $my
	
	$w.style.readout configure \
		-text [format "x: %6.2f y: %6.2f" $mx $my]
}

proc ::View3d::clear {} {
	variable w
	$w.viewport delete all
}

proc ::View3d::redraw {} {
	variable w
	
	variable points
	variable colors
	variable size
	variable zoom
	variable z0
	variable dz
	
	variable axes
	variable acolors
	variable asize
	variable azoom
	variable az0
	variable adz

	# Redraw the scene.
	$w.viewport delete all
	#$w.viewport create rectangle 0 0 $size $size -fill white
	render $points $colors [expr $size/2] [expr $size/2] \
		$zoom $z0 $dz
	render $axes $acolors [expr $size/10] [expr 9*$size/10] \
		$azoom $az0 $adz
}

proc ::View3d::drag {mx my} {
	variable w
	variable mouse_x
	variable mouse_y
	variable basis
	variable psi
	
	# Determine the angles the rotate the view.
	set phi_y [expr -$psi*($mx - $mouse_x)]
	set phi_x [expr $psi*($my - $mouse_y)]
	
	# These only commute if the angles are infinitesimal.
	set basis [transmult [transaxis y $phi_y rad] $basis]
	set basis [transmult [transaxis x $phi_x rad] $basis]
			
	redraw
	
	set mouse_x $mx
	set mouse_y $my
}



