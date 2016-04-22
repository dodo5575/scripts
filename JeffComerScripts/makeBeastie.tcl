source renderDiagram.tcl
source saveView.tcl

set step 1
set hasPdb 0
set outName sponge

if {$step == 1} {
    set sign 1
    set beastPos {-0.5 0 0}
    set lenFactor 4.2
    set brushNodes 20

    set left(A) {0.6 -0.1 -0.6}
    set left(B) {0.4 -0.2 -0.2}
    set left(C) {0.3 0.3 0.2}
    set left(D) {0.7 0.2 -0.6}

    set right(A) {0.5 -0.3 0.8}
    set right(B) {0.3 -0.2 0.1}
    set right(C) {0.2 0.3 0.6}
    set right(D) {0.8 0.1 0.6}
} elseif {$step == 2} {
    set sign 1
    set beastPos {0 0 0}
    set lenFactor 1.8
    set brushNodes 16
   
    set left(A) {1.4 0.4 -0.6}
    set left(B) {0.3 -0.8 -0.6}
    set left(C) {0.7 -0.4 -0.6}
    set left(D) {0.8 -0.7 0.1}

    set right(A) {1.1 -0.4 1.5}
    set right(B) {0.6 0.8 1.8}
    set right(C) {0.5 -1.5 -0.2}
    set right(D) {0.3 -0.3 2.6}
} elseif {$step == 3} {
    set sign -1
    set beastPos {0 0 0}
    set lenFactor 1.8
    set brushNodes 16
    
    set right(A) {1.4 0.4 -0.6}
    set right(B) {0.3 -0.8 -0.6}
    set right(C) {0.7 -0.4 -0.6}
    set right(D) {0.8 -0.7 0.1}

    set left(A) {1.1 -0.4 1.5}
    set left(B) {0.6 0.8 1.8}
    set left(C) {0.5 -1.5 -0.2}
    set left(D) {0.3 -0.3 2.6}
} elseif {$step == 4} {
    set sign -0.5
    set beastPos {1 0 0}
    set lenFactor 4.2
    set brushNodes 20
    
    set right(A) {0.6 -0.1 -0.6}
    set right(B) {0.4 -0.2 -0.2}
    set right(C) {0.3 0.3 0.2}
    set right(D) {0.7 0.2 -0.6}

    set left(A) {0.5 -0.3 0.8}
    set left(B) {0.3 -0.2 0.1}
    set left(C) {0.2 0.3 0.6}
    set left(D) {0.8 0.1 0.6}
}

set src(A) {0.06 0.5 0.0}
set src(B) {-0.05 0.7 0.7}
set src(C) {0.03 0.0 -0.5}
set src(D) {0.04 -0.7 0.7}

# Parameters:
set len 50.0
set diam 20.0
set endLen 10.0
set up {0 0 1}
set mcSteps 600

set hairN 300
set hairLen 8.0
set dx 8
set hairCount 15 
set hairRand 4.0
set ballRad [expr {0.32*$diam}]

set posListList {}
set segList {}
set beastPos [vecScale $diam $beastPos]
set z0 [expr {-0.6*$diam}]
set z1 [expr {-0.65*$diam}]

set leftX [expr {-0.5*$len}]
set rightX [expr {0.5*$len}]

foreach b {A B C D} {
    # Draw the left brush.
    foreach {x y z} $src($b) { break }
    set r0 [list [expr $leftX - $x*$diam] [expr 0.5*$y*$diam] [expr 0.5*$z*$diam]]
    foreach {x y z} $left($b) { break }
    set r1 [list [expr $leftX - $x*$len] [expr $y*$len] [expr $z*$diam]]

    set r0 [vecAdd $r0 $beastPos]
    set r1 [vecAdd $r1 $beastPos]
    
    lappend posListList [monteCarloGood [makeChain $r1 $r0 $lenFactor $brushNodes $up] $mcSteps $z1]
    lappend segList "L${b}"


    # Draw the right brush.
    foreach {x y z} $src($b) { break }
    set r0 [list [expr $rightX + $x*$diam] [expr -0.5*$y*$diam] [expr 0.5*$z*$diam]]
    foreach {x y z} $right($b) { break }
    set r1 [list [expr $rightX + $x*$len] [expr $y*$len] [expr $z*$diam]]

    set r0 [vecAdd $r0 $beastPos]
    set r1 [vecAdd $r1 $beastPos]
    
    lappend posListList [monteCarloGood [makeChain $r1 $r0 $lenFactor $brushNodes $up] $mcSteps $z1]
    lappend segList "R${b}"
}

# Make hair.
set j 0
for {set x [expr -0.6*$len]} {$x < 0.6*$len} {set x [expr {$x+$dx}]} {
    for {set i 0} {$i < $hairCount} {incr i} {
	if {$j % 2 == 0} { 
	    set theta [expr {(2.0*$pi*$i)/$hairCount}]
	} else {
	    set theta [expr {(2.0*$pi*($i+0.5))/$hairCount}]
	}
	set y [expr {cos($theta)}]
	set z [expr {sin($theta)}]

	set r0 [list $x [expr {0.5*$diam*$y}] [expr {0.5*$diam*$z}]]
	set r1 [list [expr {$x+$hairRand*rand()}] [expr {(0.5*$diam+$hairLen-$hairRand*rand())*$y}] [expr {(0.5*$diam+$hairLen-$hairRand*rand())*$z}]]
	set r0 [vecAdd $r0 $beastPos]
	set r1 [vecAdd $r1 $beastPos]
	lappend posListList [monteCarloSteps [makeChain $r1 $r0 1.1 3 $up] 5]
	lappend segList "HR${i}"
    }
    incr j
}

if {!$hasPdb} {
    makeComboTubePdb ${outName}${step}.pdb $posListList $segList
}

# Load the molecule.
mol new ${outName}${step}.pdb
source representBeastie.tcl

# Draw the walker.
setMaterial Opaque yellow
drawPill $beastPos $len $diam

# Draw the charges.
setMaterial Opaque black

if {1} {
# Draw the positive balls.
set sel [atomselect top "segname \"L.\" and resid 0"]
set posList [$sel get {x y z}]
set rot {{0 -1 0} {1 0 0} {0 0 1}}
source renderDiagram.tcl
foreach r $posList {
    drawChargeBall $r $rot $ballRad -1 iceblue black
}

# Draw the negative balls.
set sel [atomselect top "segname \"R.\" and resid 0"]
set posList [$sel get {x y z}]
set rot {{0 -1 0} {1 0 0} {0 0 1}}
source renderDiagram.tcl
foreach r $posList {
    drawChargeBall $r $rot $ballRad 1 pink black
}

# Add charge.
set metalZ [expr {-25.0-$ballRad}]
set x0 [expr {-3.5*$len}]
set x1 [expr {3.5*$len}]
set dx [expr {2.0*$ballRad}]
for {set x $x0} {$x < $x1} {set x [expr {$x+$dx}]} {
    set r [list $x 0 $metalZ]
    drawChargeBallClear $r $rot $ballRad $sign black
}
}
# Load the view.
loadView sponge
