# computeForceOnProtein
#
# author: dbwells2@uiuc.edu

# Procs which NAMD has defined
proc atomid {segname resid atomname} {
    set sel [atomselect top "segname $segname and resid $resid and name $atomname"]
    set index [$sel get index]
    $sel delete
    return $index
}

proc addatom {atomid} {
    return
}

proc print {text} {
    puts $text
}

proc readGeometryFile { geometryFile } {
    set inCh1		[open $geometryFile r]
    set i		0
    set halfSpacing	1.5
    set seltext ""
    foreach line [split [read $inCh1] \n] {
	if {[llength $line] == 3} {
	    foreach {z R var} $line { break }

	    if {$i == 0} {
		set zDown [expr $z - $halfSpacing]
		set seltext "DNA and ((z > [expr $z -$halfSpacing] and z < [expr $z + $halfSpacing] and x*x +y*y < $R*$R)"
	    } else {

		set seltext "$seltext or (z > [expr $z -$halfSpacing] and z < [expr $z + $halfSpacing] and x*x +y*y < $R*$R)"
	    }
	    incr i	
	}
    }
    set seltext "$seltext )"
    
    return $seltext
}

proc computeForceOnDNA { pdb psf dcd region seltext gridForceScript scale } {
    puts "Computing force on DNA using following parameters:"
    puts "	pdb		= $pdb"
    puts "	psf		= $psf"
    puts "	dcd		= $dcd"
    puts "	region		= $region"
    puts "	seltext		= $seltext"
    puts "	gridForceScript	= $gridForceScript"
    puts "	scale		= $scale"
    
    global forceScale
    global gridPotential
    global e1x e1y e1z e2x e2y e2z e3x e3y e3z
    global originX originY originZ
    global nX nY nZ
    global inv1 inv2 inv3
    
    set tramol [mol load psf $psf]
    #mol addfile $dcd last 0
    #set trasel [atomselect $tramol $seltext]
    
    #set nframes [molinfo top get numframes]
    set nframes [exec catdcd -num $dcd | grep Total | awk "{print \$3}"]
    
    set outCh [open ${dcd}.${region}.DNAForce "w"]
    
    # Gridforce parameters
    set targetAtomPdb		$pdb
    set gridPotentialFile	/home/david/Work/polyA-gf-IV/tmp.dx
    set forceScale		"0.0 0.0 ${scale}"
    set forcesRecalcFreq	1	;# doesn't actually affect anything here
    
    source $gridForceScript
    
    for { set frame 0 } { $frame < $nframes } { incr frame } {
	mol addfile $dcd molid $tramol first $frame last $frame
	set trasel [atomselect $tramol $seltext]
	
	set trapos [$trasel get {x y z}]
	set charges [$trasel get charge]
	
	set totalForce {0 0 0}
	foreach tpos $trapos q $charges {
	    foreach {x y z} $tpos { break }
	    set force_kcal [getGridForce $x $y $z $q]
	    set force [vecscale $force_kcal 69.523]
	    
	    puts "force = $force"
	    
	    set totalForce [vecadd $totalForce $force]
	}
	
	puts $outCh "$frame $totalForce"
	puts [format "FRAME %5s/%-5s  force = %s pN" $frame [expr { $nframes - 1 }] $totalForce]
	
	animate delete all
	$trasel delete
    }
    
    close $outCh
}
