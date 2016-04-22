#
# This script read a DCD with many frames and write a file with the connectivity
# with the beta values change
# useful for movies
#
#  vmd -dispdev win -e WRAPDCD.tcl




##################### INPUT - OUTPUT VALUES ################################

set psfFile sio2crystal.psf
set dcdFile simALL.dcd
set output  simALL.wrap

set DCDcenter "27.38  27.38  27.50"

# Cell vectors

set pbcX  "54.758 0 0"
set pbcY  "0 54.758 0"
set pbcZ "0 0 55.584"



########################## MAIN PART STARTS HERE ##########################




#  General Info
# ==============


# Load molecule
mol load psf $psfFile dcd $dcdFile

# Number of steps
set num_steps [molinfo top get numframes]


# Calculate number of cell to replicate

#### INPUT ####





#set diag [ vecadd $pbcX $pbcY $pbcZ ]

proc CellImages { initialPoint diag min max } {


    # find border points of the periodic cell

    foreach { XminCell YminCell ZminCell } $initialPoint { break }
    foreach { XLength YLength ZLength }  $diag { break }

    set  XmaxCell [ expr  $XLength + $XminCell]
    set  YmaxCell [ expr  $YLength + $YminCell]
    set  ZmaxCell [ expr  $ZLength + $ZminCell]


    # find border points of the frame

    foreach { XminFrame YminFrame ZminFrame } $min { break }
    foreach { XmaxFrame YmaxFrame ZmaxFrame } $max { break }


  
    # how many cell in the positive directions

    set tmpx $XmaxCell
    set tmpy $YmaxCell
    set tmpz $ZmaxCell

    set i 0
 
    while { $tmpx  < $XmaxFrame  } {

	incr i

	set borderLeft $tmpx
	set moveX [ expr $XLength * $i * (-1) ]
	set tmpx  [ expr $tmpx + $XLength ]
	set borderRight $tmpx

	lappend arrayX "$borderLeft $borderRight $moveX"

    }

    set i 0

    while { $tmpy  < $YmaxFrame  } {

	incr i

	set borderLeft $tmpy
	set moveY [ expr $YLength * $i * (-1)  ]
	set tmpy  [ expr $tmpy + $YLength ]
	set borderRight $tmpy

	lappend arrayY "$borderLeft $borderRight $moveY"

    }

    set i 0

    while { $tmpz  < $ZmaxFrame  } {

	incr i

	set borderLeft $tmpz
	set moveZ [ expr $ZLength * $i  * (-1) ]
	set tmpz  [ expr $tmpz + $ZLength ]
	set borderRight $tmpz

	lappend arrayZ "$borderLeft $borderRight $moveZ"

    }



    # how many cell in the negative directions

    set tmpx $XminCell
    set tmpy $YminCell
    set tmpz $ZminCell

    set i 0
 
    while { $tmpx  > $XminFrame  } {

	incr i

	set borderRight $tmpx
	set moveX [ expr $XLength * $i ]
	set tmpx  [ expr $tmpx - $XLength ]
	set borderLeft $tmpx

	lappend arrayX "$borderLeft $borderRight $moveX"

    }

    set i 0

    while { $tmpy  > $YminFrame  } {

	incr i

	set borderRight $tmpy
	set moveY [ expr $YLength * $i  ]
	set tmpy  [ expr $tmpy - $YLength ]
	set borderLeft $tmpy

	lappend arrayY "$borderLeft $borderRight $moveY"

    }

    set i 0

    while { $tmpz  > $ZminFrame  } {

	incr i

	set borderRight $tmpz
	set moveZ [ expr $ZLength * $i  ]
	set tmpz  [ expr $tmpz - $ZLength ]
	set borderLeft $tmpz

	lappend arrayZ "$borderLeft $borderRight $moveZ"

    }


    if ![ info exist arrayX ] {
	set arrayX none
    } 

    if ![ info exist arrayY ] {
	set arrayY none
    }


    if ![ info exist arrayZ ] {
	set arrayZ none
    }


    lappend arrayALL $arrayX $arrayY $arrayZ

    return $arrayALL

}


#################


set k 0

for {set frame 0} {$frame < $num_steps} {incr frame} {

    puts "wrapping frame $k"

    ### 1 Save working frame

    #Select frame
    set currentFrame [atomselect top "all" frame $k]

    
    # Measure MinMax of frame
    foreach {min max} [ measure minmax $currentFrame ] { break }

    # Measure PBC points
    set diag [ vecadd $pbcX $pbcY $pbcZ ]
    set 2center [ vecscale  2  $DCDcenter ]
    set initialPoint [ vecsub $2center $diag ] 

 
    foreach { a b c }  [ CellImages $initialPoint $diag $min $max ] { break }

  
    if { $a != "none" } {

	foreach xelement $a {
	    foreach { Left Right Move } $xelement { break }
	    set outofframe [atomselect top "x > $Left and x <= $Right" frame $k]
	    $outofframe moveby "$Move 0 0"
	    $outofframe delete    
	}
	
    } 
    
    if { $b != "none" } {

	foreach yelement $b {
	    foreach { Left Right Move } $yelement { break }
	    set outofframe [atomselect top "y > $Left and y <= $Right" frame $k]
	    $outofframe moveby "0 $Move 0"
	    $outofframe delete
	}

    }

    if { $c != "none" } {

	foreach zelement $c {
	    foreach { Left Right Move } $zelement { break }
	    set outofframe [atomselect top "z > $Left and z <= $Right" frame $k]
	    $outofframe moveby "0 0 $Move"
	    $outofframe delete    
	}
	
    }


    $currentFrame writepdb temp($k).pdb
    $currentFrame delete 
    incr k
}
 
mol delete top

### Write DCD 

puts "writing dcd"

set i 0

mol load psf $psfFile
set MolID [molinfo top]

while { $i < $num_steps } {

    mol addfile temp($i).pdb type {pdb} first 0 last -1 step 1 waitfor 1  $MolID
    animate style Loop
    exec rm temp($i).pdb
    incr i
	
}

set nf [ expr $num_steps - 1 ]

set sel [atomselect top "all"]
animate write dcd $output.dcd beg 0 end $nf waitfor all sel $sel skip 1 $MolID
$sel delete

mol delete top

### Visualize DCD 

#load molecule
puts "loading molecule"
mol load psf $psfFile dcd $output.dcd

#draw periodic cell

foreach { minx miny minz  }  $initialPoint { break } 

foreach { maxx maxy maxz  }  [ vecadd $diag $initialPoint ]  { break }


 
# and draw the lines
draw materials off
draw color blue
draw line "$minx $miny $minz" "$maxx $miny $minz" width 5
draw line "$minx $miny $minz" "$minx $maxy $minz" width 5
draw line "$minx $miny $minz" "$minx $miny $maxz" width 5

draw line "$maxx $miny $minz" "$maxx $maxy $minz" width 5
draw line "$maxx $miny $minz" "$maxx $miny $maxz" width 5

draw line "$minx $maxy $minz" "$maxx $maxy $minz" width 5
draw line "$minx $maxy $minz" "$minx $maxy $maxz" width 5

draw line "$minx $miny $maxz" "$maxx $miny $maxz" width 5
draw line "$minx $miny $maxz" "$minx $maxy $maxz" width 5

draw line "$maxx $maxy $maxz" "$maxx $maxy $minz" width 5
draw line "$maxx $maxy $maxz" "$minx $maxy $maxz" width 5
draw line "$maxx $maxy $maxz" "$maxx $miny $maxz" width 5




