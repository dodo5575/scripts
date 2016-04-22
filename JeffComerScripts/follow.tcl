# Author: Jeff Comer <jcomer2@illinois.edu>
# Parameters:
set prefix frames/tachyon_40ps
set keyInterval 20
set stride 1
set render 1
set selText "segname A60 and resid 0"
set followPos {13 0 0}
set dcdFreq 5000
set lx 113.394
set ly 113.402
set lz 56.6647
set subjectMol 3

#set dt [expr 0.1*5000/$dcdFreq]
set dt 0.1
display projection perspective

proc setView0 {origin d} {
    if {[vecdot $d $d] == 0.0} { set d [list 0 0 1] }
    set up [list 0.1 0.9 0.0]
    set up [vecscale [expr 1.0/[veclength $up]] $up]
    set ez [vecscale [expr -1.0/[veclength $d]] $d]
    set ey [vecsub $up [vecscale [vecdot $up $ez] $ez]]
    set ey [vecscale [expr 1.0/[veclength $ey]] $ey]
    set ex [veccross $ey $ez]

    set row0 [list [lindex $ex 0] [lindex $ey 0] [lindex $ez 0] 0.0]
    set row1 [list [lindex $ex 1] [lindex $ey 1] [lindex $ez 1] 0.0]
    set row2 [list [lindex $ex 2] [lindex $ey 2] [lindex $ez 2] 0.0]
    set row3 [list 0.0 0.0 0.0 1.0]
    set basis [list $row0 $row1 $row2 $row3]
    #set roll [list {0 -1 0 0} {1 0 0 0} {0 0 1 0} {0 0 0 1}]
    set roll [transaxis z -90]
    set basis [transmult $roll $basis]
}

proc setView {origin d} {
    display distance -2.0
    display focallength 2.0
    set s 0.075

    #set dir [list [lindex $d 0] [lindex $d 1] 0]
    #set dir [vecscale [expr 1.0/[veclength $dir] $dir]

    set rotMat [list {0 1 0 0} {0 0 1 0} {1 0 0 0} {0 0 0 1}]
    #set rotMat $basis

    set ox [expr -[lindex $origin 0] + 2.0/$s]
    set oy [expr -[lindex $origin 1]]
    set oz [expr -[lindex $origin 2]]
    set cenMat [list [list 1 0 0 $ox] [list 0 1 0 $oy] [list 0 0 1 $oz] [list 0 0 0 1]]

    set scaleMat [list [list $s 0 0 0] [list 0 $s 0 0] [list 0 0 $s 0] [list 0 0 0 1]]
    set globeMat [transidentity]

    set viewMat [list $cenMat $rotMat $scaleMat $globeMat]
    set molList [molinfo list]
    foreach m $molList {
	molinfo $m set {center_matrix rotate_matrix scale_matrix global_matrix} $viewMat
    }
}

proc renderFrameSnap {prefix} {
    render snapshot $prefix.tga ""
}

proc renderFrameTachyon {prefix} {
    render Tachyon $prefix.dat ""
}

proc wrap {d} {
    global lx ly lz

    foreach {dx dy dz} $d { break }
    while {$dx > [expr 0.5*$lx]} { set dx [expr $dx - $lx] }
    while {$dx < [expr -0.5*$lx]} { set dx [expr $dx + $lx] }
    while {$dy > [expr 0.5*$ly]} { set dy [expr $dy - $ly] }
    while {$dy < [expr -0.5*$ly]} { set dy [expr $dy + $ly] }
    while {$dz > [expr 0.5*$lz]} { set dz [expr $dz - $lz] }
    while {$dz < [expr -0.5*$lz]} { set dz [expr $dz + $lz] }
    
    return [list $dx $dy $dz]
}

proc wrap1 {d} {
    global lx ly lz

    foreach {dx dy dz} $d { break }
    while {$dx > [expr 0.5*$lx]} { set dx [expr $dx - $lx] }
    while {$dx < [expr -0.5*$lx]} { set dx [expr $dx + $lx] }
    while {$dy > [expr 0.5*$ly]} { set dy [expr $dy - $ly] }
    while {$dy < [expr -0.5*$ly]} { set dy [expr $dy + $ly] }
    #while {$dz > [expr 0.5*$lz]} { set dz [expr $dz - $lz] }
    #while {$dz < [expr -0.5*$lz]} { set dz [expr $dz + $lz] }
    while {$dz > $lz} { set dz [expr $dz - $lz] }
    while {$dz < 0.0} { set dz [expr $dz + $lz] }
    
    return [list $dx $dy $dz]
}

proc goalPosition1 {pos vel followDist r} {
    global lx ly lz
    
    if {[vecdot $vel $vel] == 0.0} { set vel [list 0 0 1] }
    set v [vecscale [expr $followDist/[veclength $vel]] $vel]

    set d [vecsub [vecsub $pos $v] $r]
    return [vecadd $r [wrap $d]]
}

proc goalPosition {pos vel followPos r} {
    global lx ly lz

    set d [vecsub [vecsub $pos $followPos] $r]
    return [vecadd $r [wrap $d]]
}

proc drawSphere {pos} {
    graphics top sphere $pos radius 1.0 resolution 25
    return
}

proc insertControlPoints { vert } {
    set last [expr [llength $vert]-1]
    set ret {}
    
    # Set the control point spacing parameter.
    set alpha [expr 1.0/3.0]
    
    # Determine the control points neighboring the internal vertices.		
    for {set i 1} {$i < $last} {incr i} {
	# Determine the control points.
	set ra [lindex $vert [expr $i-1]]
	set r [lindex $vert $i]
	set rb [lindex $vert [expr $i+1]]
	
	set amag [veclength [vecsub $ra $r]]
	set bmag [veclength [vecsub $rb $r]]		
	set aunit [vecnorm [vecsub $ra $r]]
	set bunit [vecnorm [vecsub $rb $r]]
	set side [vecnorm [vecsub $bunit $aunit]]
	
	# These are the control points that neighbor vertex[i].
	set ctrlA [vecadd $r [vecscale $side [expr -$alpha*$amag]]]
	set ctrlB [vecadd $r [vecscale $side [expr $alpha*$bmag]]]
	
	# Add the new points to "ret".
	lappend ret $ctrlA $r $r $ctrlB
    }
    
    # Set the first control point linearly.
    set r [lindex $vert 0]
    set dmag [veclength [vecsub [lindex $ret 0] $r]] 
    set side [vecnorm [vecsub [lindex $ret 0] $r]]
    set ctrlB [vecadd $r [vecscale $side [expr $alpha*$dmag]]]
    set ret [concat [list $r] [list $ctrlB] $ret]
    
    # Set the last control point linearly.
    set r [lindex $vert $last]
    set d [vecsub $r [lindex $ret [expr [llength $ret]-1]]]
    set dmag [veclength $d]
    set side [vecnorm $d]
    set ctrlA [vecadd $r [vecscale $side [expr -$alpha*$dmag]]]
    lappend ret $ctrlA $r
    
    # Return the vertex/control point list.
    return $ret
}

# Compute the Bezier spline vectors for each vertex.
# Return sequences of four vectors a, b, c and p0. 
proc bezierParameters { point } {
    set ret {}
    
    foreach {p0 p1 p2 p3} $point {
	set c [vecscale [vecsub $p1 $p0] 3.0]
	set b [vecsub [vecscale [vecsub $p2 $p1] 3.0] $c]
	set a [vecadd $p3 [vecinvert $p0] [vecinvert $c] [vecinvert $b]]
	
	lappend ret $a $b $c $p0
    }
    return $ret
}

# Generate points on a spline given Bezier parameters of vertices.
proc sampleSpline { samples param } {
    set ret {}
    
    foreach {a b c p0} $param {
	for {set i 0} {$i < $samples} {incr i} {
	    set t [expr 1.0*$i/$samples]
	    set p1 [vecscale $c $t]
	    set p2 [vecscale $b [expr $t*$t]]
	    set p3 [vecscale $a [expr $t*$t*$t]]
	    
	    # Locate the point along the spline.
	    set r [vecadd $p0 $p1 $p2 $p3]
	    
	    lappend ret $r
	}
    }
    return $ret	
}

# Get the unit tangent vector for a curve.
proc unitTangent { curve } {
    set ret {}

    # Make the first and last vectors.
    set v0 [vecsub [lindex $curve 1] [lindex $curve 0]]
    set v0 [vecscale [expr 1.0/[veclength $v0]] $v0]
    set v1 [vecsub [lindex $curve end] [lindex $curve end-1]]
    set v1 [vecscale [expr 1.0/[veclength $v1]] $v1]

    # Make the unit tangent vectors.
    lappend ret $v0
    set n1 [expr [llength $curve] - 1]
    for {set i 1} {$i < $n1} {incr i} {
	set i0 [expr $i-1]
	set i1 [expr $i+1]
	set v [vecsub [lindex $curve $i0] [lindex $curve $i1]]
	set v [vecscale [expr 1.0/[veclength $v]] $v]

	lappend ret $v
    }
    lappend ret $v1

    return $ret
}

#####################################
## Main portion of script
set nFrames [molinfo top get numframes]
set sel [atomselect top $selText]
set subject [atomselect $subjectMol all]
graphics top delete all

# Get the frame positions, unwrapping.
set framePos {}
molinfo top set frame 0
set pos0 [measure center $sel weight mass]
lappend framePos $pos0
for {set f 1} {$f < $nFrames} {incr f} {
    molinfo top set frame $f
    set pos [measure center $sel weight mass]
    set dr [wrap [vecsub $pos $pos0]]
    set pos1 [vecadd $pos0 $dr]

    lappend framePos $pos1
    set pos0 $pos1
}

# Compute the key frame positions.
set path {}
for {set f 0} {$f < $nFrames} {incr f $keyInterval} {
    lappend path [lindex $framePos $f]
}
lappend path [lindex $framePos end]

# Make a spline connecting the key frames.
set param [bezierParameters [insertControlPoints $path]]
set spline [sampleSpline $keyInterval $param]
set tangent [unitTangent $spline]

# Draw the spline.
#for {set f 0} {$f < $nFrames} {incr f} {
#    molinfo top set frame $f
#    graphics top color 0
#    drawSphere [lindex $framePos $f]
#    graphics top color 1
    #drawSphere [lindex $spline $f]
#    drawSphere [wrap1 [vecadd [lindex $spline $f] $followPos]]

#    display update
#}

for {set f 0} {$f < $nFrames} {incr f $stride} {
    molinfo top set frame $f

    set splinePos [wrap [lindex $spline $f]]
    
    # Move the molecule to this position.
    #$subject moveby [vecinvert [measure center $subject weight mass]]
    #$subject moveby $splinePos

    set camPos [wrap1 [vecadd [lindex $spline $f] $followPos]]
    #set camVel [lindex $tangent $f]
    set camVel [wrap [vecsub $camPos [lindex $framePos $f]]]

    setView $camPos $camVel

    display update
    if {$render} {
	renderFrameTachyon ${prefix}${f}
    }
    puts "frame $f"    
}
