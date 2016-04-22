source renderDiagram.tcl

set step 4
set outName beastie

if {$step == 1} {
    set sign 1
    set beastPos {-1 0 0}
    set lenFactor 4
    set brushNodes 12

    set src(A) {0.1 0.1 0.2}
    set src(B) {-0.05 0.1 0.1}
    set src(C) {0.1 -0.2 0.3}

    set left(A) {0.6 -0.1 -0.6}
    set left(B) {0.4 -0.2 -0.6}
    set left(C) {0.3 0.3 -0.6}

    set right(A) {0.5 -0.1 0.8}
    set right(B) {0.6 0.2 0.4}
    set right(C) {0.2 0.3 0.6}
} elseif {$step == 2} {
    set sign 1
    set beastPos {0 0 0}
    set lenFactor 1.1
    set brushNodes 5

    set src(A) {0.1 0.1 0.2}
    set src(B) {-0.05 0.1 0.1}
    set src(C) {0.1 -0.2 0.3}

    set left(A) {1.0 -0.4 -0.6}
    set left(B) {0.2 -0.8 -0.6}
    set left(C) {0.6 0.4 -0.6}

    set right(A) {0.9 -0.4 1.5}
    set right(B) {0.3 1.1 2.0}
    set right(C) {0.6 -1.0 1.0}
} elseif {$step == 3} {
    set sign -1
    set beastPos {0 0 0}
    set lenFactor 1.1
    set brushNodes 5

    set src(A) {0.1 0.1 0.2}
    set src(B) {-0.05 0.1 0.1}
    set src(C) {0.1 -0.2 0.3}

    set right(A) {1.0 -0.4 -0.6}
    set right(B) {0.2 -0.8 -0.6}
    set right(C) {0.6 0.4 -0.6}

    set left(A) {0.9 -0.4 1.5}
    set left(B) {0.3 1.1 2.0}
    set left(C) {0.6 -1.0 1.0}
} elseif {$step == 4} {
    set sign -1
    set beastPos {1 0 0}
    set lenFactor 4
    set brushNodes 12

    set src(A) {0.1 0.1 0.2}
    set src(B) {-0.05 0.1 0.1}
    set src(C) {0.1 -0.2 0.3}

    set right(A) {0.6 -0.1 -0.6}
    set right(B) {0.4 -0.2 -0.6}
    set right(C) {0.3 0.3 -0.6}

    set left(A) {0.5 -0.1 0.8}
    set left(B) {0.6 0.2 0.4}
    set left(C) {0.2 0.3 0.6}
}

set len 50.0
set diam 20.0
set endLen 10.0
set mcSteps 400
set up {0 0 1}

set posListList {}
set segList {}
set hairN 300
set hairLen 8.0

set beastPos [vecScale $diam $beastPos]
set z0 [expr {-0.6*$diam}]
set z1 [expr {-0.65*$diam}]

set leftX [expr {-0.5*$len}]
set rightX [expr {0.5*$len}]

# Make hair.
set dx 8
set hairCount 15 
set hairRand 4.0
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

# Load the molecule.
mol new ${outName}${step}.pdb
source representBeastie.tcl

# Draw the walker.
setMaterial Opaque yellow
drawPill $beastPos $len $diam

setMaterial Opaque black

# Draw the positive balls.
set sel [atomselect top "segname LA LB LC and resid 0"]
set posList [$sel get {x y z}]
set rot {{0 -1 0} {1 0 0} {0 0 1}}
source renderDiagram.tcl
foreach r $posList {
    drawChargeBall $r $rot [expr {0.25*$diam}] -1 iceblue black
}

# Draw the negative balls.
set sel [atomselect top "segname RA RB RC and resid 0"]
set posList [$sel get {x y z}]
set rot {{0 -1 0} {1 0 0} {0 0 1}}
source renderDiagram.tcl
foreach r $posList {
    drawChargeBall $r $rot [expr {0.25*$diam}] 1 iceblue black
}

# Add charge.
set ballRad [expr {0.25*$diam}]
set metalZ [expr {-25.0-$ballRad}]
set x0 [expr {-3.5*$len}]
set x1 [expr {3.5*$len}]
set dx [expr {2.0*$ballRad}]
for {set x $x0} {$x < $x1} {set x [expr {$x+$dx}]} {
    set r [list $x 0 $metalZ]
    drawChargeBallClear $r $rot $ballRad $sign black
}

# Load the view.
loadView pill
