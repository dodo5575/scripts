##########################################################
##This script produces different view of the trojectory.
## 
##Please read the "Readme" file for description of the commands in this script.
##
##Command List:
##    view1
##    view2
##    view3
##    view4
##    viewAT
##    viewCG
##    rotateX
##    rotateY
##    rotateX
##    zoomIn
##    zoomOut
##    Play
##    Pause
##    Rewind
##    animation
##
##written by Chen-Yu Li
##cli56@illinois.edu
##2013/7/18
##########################################################


proc settings {} {
	color Display Background white
	display projection Orthographic
	display nearclip set 0.01
	axes location Off
}


proc reSet {} {

	display resetview
	settings

	set NumOfRep [lindex [mol list top] 12]
	
	for {set i 1} {$i <= $NumOfRep} {incr i} {
		mol delrep 0 top
	}

}


proc getTopMolID {} {

	set id [molinfo top]
	return $id
}


proc view1 {} {

	reSet
	reSet
	
	set id [getTopMolID]

	mol addrep $id
	mol modcolor 0 $id ColorID 21
	mol modstyle 0 $id VDW 0.400000 12.000000	
	mol modselect 0 $id nucleic and backbone and name P
	mol material Opaque

	mol addrep $id
	mol modcolor 1 $id SegName
	mol modstyle 1 $id DynamicBonds 1.600000 0.300000 6.000000
	mol modselect 1 $id nucleic and not backbone
	mol material Opaque
	
	#scale by 0.5

}


proc view2 {} {

	reSet
	reSet

	set id [getTopMolID]
	
	mol addrep $id
	mol modcolor 0 $id ColorID 21
	mol modstyle 0 $id Polyhedra 1.6
	mol modselect 0 $id nucleic and backbone and (not segname SC0 SC1 SC2 SC3)
	mol material Opaque

	mol addrep $id
	mol modcolor 1 $id ColorID 1
	mol modstyle 1 $id Polyhedra 1.6
	mol modselect 1 $id nucleic and backbone and segname SC0 SC1 SC2 SC3
	mol material Opaque

	mol addrep $id
	mol modcolor 2 $id ResName
	mol modstyle 2 $id CPK 1.000000 0.300000 10.000000 10.000000
	mol modselect 2 $id nucleic and not backbone
	mol material Opaque

	#scale by 0.5

}

proc view3 {} {

	reSet
	reSet

	set id [getTopMolID]
	
	mol addrep $id
	mol modcolor 0 $id ColorID 21
	mol modstyle 0 $id QuickSurf 0.700000 0.970000 1.940000 0.000000
	mol modselect 0 $id nucleic and backbone
	mol modmaterial 0 $id Transparent

	mol addrep $id
	mol modcolor 1 $id ColorID 13
	mol modstyle 1 $id Bonds 0.300000 10.000000
	mol modselect 1 $id nucleic and not backbone
	mol modmaterial 1 $id AOChalky

	mol addrep $id
	mol modcolor 2 $id Name
	mol modstyle 2 $id Licorice 0.200000 50.000000 50.000000
	mol modselect 2 $id segname MGHH
	mol modmaterial 2 $id Opaque

	mol addrep $id
	mol modcolor 3 $id Name
	mol modstyle 3 $id VDW 0.600000 12.000000
	mol modselect 3 $id segname MGHH and name MG
	mol modmaterial 3 $id Opaque

	mol addrep $id
	mol modcolor 4 $id Name
	mol modstyle 4 $id Lines 1.000000
	mol modselect 4 $id water and noh
	mol modmaterial 4 $id Opaque
	
	#scale by 0.5
	rotate x by 90

}


proc view4 {} {

	reSet
	reSet

	set id [getTopMolID]
	
	mol addrep $id
	mol modcolor 0 $id ColorID 23
	mol modstyle 0 $id VDW 0.5 12
	mol modselect 0 $id nucleic and backbone
	mol modmaterial 0 $id Opaque

	mol addrep $id
	mol modcolor 1 $id ColorID 3
	mol modstyle 1 $id Licorice 0.200000 50.000000 50.000000
	mol modselect 1 $id nucleic and not backbone
	mol modmaterial 1 $id AOChalky

	mol addrep $id
	mol modcolor 2 $id Name
	mol modstyle 2 $id VDW 0.5 12
	mol modselect 2 $id ion or name MG
	mol modmaterial 2 $id Glossy

	mol addrep $id
	mol modcolor 3 $id ColorID 0
	mol modstyle 3 $id Lines 1.000000
	mol modselect 3 $id water and noh
	mol modmaterial 3 $id Opaque
	
	#scale by 0.5
	rotate x by 90

}

proc viewCG {} {
	reSet
	reSet

	set id [getTopMolID]

	set viewplist {}
	set fixedlist {}
	mol representation Licorice 0.200000 50.000000 50.000000
	mol color Name
	mol selection {(segname P1 and resid 30) or (segname SC1 and resid 15)}
	mol material Opaque
	mol addrep $id
	mol selupdate 0 $id 0
	mol colupdate 0 $id 0
	mol scaleminmax $id 0 0.000000 0.000000
	mol smoothrep $id 0 0
	mol drawframes $id 0 {now}
	mol clipplane center 0 0 $id {0.0 0.0 0.0}
	mol clipplane color  0 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 0 0 $id {0.0 0.0 1.0}
	mol clipplane status 0 0 $id {0}
	mol clipplane center 1 0 $id {0.0 0.0 0.0}
	mol clipplane color  1 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 1 0 $id {0.0 0.0 1.0}
	mol clipplane status 1 0 $id {0}
	mol clipplane center 2 0 $id {0.0 0.0 0.0}
	mol clipplane color  2 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 2 0 $id {0.0 0.0 1.0}
	mol clipplane status 2 0 $id {0}
	mol clipplane center 3 0 $id {0.0 0.0 0.0}
	mol clipplane color  3 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 3 0 $id {0.0 0.0 1.0}
	mol clipplane status 3 0 $id {0}
	mol clipplane center 4 0 $id {0.0 0.0 0.0}
	mol clipplane color  4 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 4 0 $id {0.0 0.0 1.0}
	mol clipplane status 4 0 $id {0}
	mol clipplane center 5 0 $id {0.0 0.0 0.0}
	mol clipplane color  5 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 5 0 $id {0.0 0.0 1.0}
	mol clipplane status 5 0 $id {0}
	mol representation HBonds 3.500000 25.000000 5.000000
	mol color Name
	mol selection {(segname P1 and resid 30) or (segname SC1 and resid 15)}
	mol material Opaque
	mol addrep $id
	mol selupdate 1 $id 0
	mol colupdate 1 $id 0
	mol scaleminmax $id 1 0.000000 0.000000
	mol smoothrep $id 1 0
	mol drawframes $id 1 {now}
	mol clipplane center 0 1 $id {0.0 0.0 0.0}
	mol clipplane color  0 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 0 1 $id {0.0 0.0 1.0}
	mol clipplane status 0 1 $id {0}
	mol clipplane center 1 1 $id {0.0 0.0 0.0}
	mol clipplane color  1 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 1 1 $id {0.0 0.0 1.0}
	mol clipplane status 1 1 $id {0}
	mol clipplane center 2 1 $id {0.0 0.0 0.0}
	mol clipplane color  2 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 2 1 $id {0.0 0.0 1.0}
	mol clipplane status 2 1 $id {0}
	mol clipplane center 3 1 $id {0.0 0.0 0.0}
	mol clipplane color  3 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 3 1 $id {0.0 0.0 1.0}
	mol clipplane status 3 1 $id {0}
	mol clipplane center 4 1 $id {0.0 0.0 0.0}
	mol clipplane color  4 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 4 1 $id {0.0 0.0 1.0}
	mol clipplane status 4 1 $id {0}
	mol clipplane center 5 1 $id {0.0 0.0 0.0}
	mol clipplane color  5 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 5 1 $id {0.0 0.0 1.0}
	mol clipplane status 5 1 $id {0}
	#mol rename $id square2plate-1MKCl.psf
	set viewpoints([molinfo top]) {{{1 0 0 1.12573} {0 1 0 0.236681} {0 0 1 0.164417} {0 0 0 1}} {{-0.00595021 -0.903083 0.429456 0} {0.293062 0.40903 0.864196 0} {-0.956081 0.130999 0.262224 0} {0 0 0 1}} {{0.0330961 0 0 0} {0 0.0330961 0 0} {0 0 0.0330961 0} {0 0 0 1}} {{1 0 0 0.45} {0 1 0 0.22} {0 0 1 0} {0 0 0 1}}}
	lappend viewplist [molinfo top]
	set topmol [molinfo top]
	# done with molecule 0
	foreach v $viewplist {
	  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)
	}
	foreach v $fixedlist {
	  molinfo $v set fixed 1
	}
	unset viewplist
	unset fixedlist
	#mol $id $topmol
	unset topmol
	
}

proc viewAT {} {

	reSet
	reSet

	set id [getTopMolID]

	set viewplist {}
	set fixedlist {}
	mol representation Licorice 0.200000 50.000000 50.000000
	mol color Name
	mol selection {(segname SC1 and resid 21) or (segname P1 and resid 24)}
	mol material Opaque
	mol addrep top
	mol selupdate 0 $id 0
	mol colupdate 0 $id 0
	mol scaleminmax $id 0 0.000000 0.000000
	mol smoothrep $id 0 0
	mol drawframes $id 0 {now}
	mol clipplane center 0 0 $id {0.0 0.0 0.0}
	mol clipplane color  0 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 0 0 $id {0.0 0.0 1.0}
	mol clipplane status 0 0 $id {0}
	mol clipplane center 1 0 $id {0.0 0.0 0.0}
	mol clipplane color  1 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 1 0 $id {0.0 0.0 1.0}
	mol clipplane status 1 0 $id {0}
	mol clipplane center 2 0 $id {0.0 0.0 0.0}
	mol clipplane color  2 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 2 0 $id {0.0 0.0 1.0}
	mol clipplane status 2 0 $id {0}
	mol clipplane center 3 0 $id {0.0 0.0 0.0}
	mol clipplane color  3 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 3 0 $id {0.0 0.0 1.0}
	mol clipplane status 3 0 $id {0}
	mol clipplane center 4 0 $id {0.0 0.0 0.0}
	mol clipplane color  4 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 4 0 $id {0.0 0.0 1.0}
	mol clipplane status 4 0 $id {0}
	mol clipplane center 5 0 $id {0.0 0.0 0.0}
	mol clipplane color  5 0 $id {0.5 0.5 0.5 }
	mol clipplane normal 5 0 $id {0.0 0.0 1.0}
	mol clipplane status 5 0 $id {0}
	mol representation HBonds 3.500000 24.000000 5.000000
	mol color Name
	mol selection {(segname SC1 and resid 21) or (segname P1 and resid 24)}
	mol material Opaque
	mol addrep $id
	mol selupdate 1 $id 0
	mol colupdate 1 $id 0
	mol scaleminmax $id 1 0.000000 0.000000
	mol smoothrep $id 1 0
	mol drawframes $id 1 {now}
	mol clipplane center 0 1 $id {0.0 0.0 0.0}
	mol clipplane color  0 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 0 1 $id {0.0 0.0 1.0}
	mol clipplane status 0 1 $id {0}
	mol clipplane center 1 1 $id {0.0 0.0 0.0}
	mol clipplane color  1 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 1 1 $id {0.0 0.0 1.0}
	mol clipplane status 1 1 $id {0}
	mol clipplane center 2 1 $id {0.0 0.0 0.0}
	mol clipplane color  2 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 2 1 $id {0.0 0.0 1.0}
	mol clipplane status 2 1 $id {0}
	mol clipplane center 3 1 $id {0.0 0.0 0.0}
	mol clipplane color  3 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 3 1 $id {0.0 0.0 1.0}
	mol clipplane status 3 1 $id {0}
	mol clipplane center 4 1 $id {0.0 0.0 0.0}
	mol clipplane color  4 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 4 1 $id {0.0 0.0 1.0}
	mol clipplane status 4 1 $id {0}
	mol clipplane center 5 1 $id {0.0 0.0 0.0}
	mol clipplane color  5 1 $id {0.5 0.5 0.5 }
	mol clipplane normal 5 1 $id {0.0 0.0 1.0}
	mol clipplane status 5 1 $id {0}
	#mol rename top square2plate-1MKCl.psf
	set viewpoints([molinfo top]) {{{1 0 0 1.12573} {0 1 0 0.236681} {0 0 1 0.164417} {0 0 0 1}} {{0.063384 0.920929 -0.384577 0} {0.237702 -0.388189 -0.89041 0} {-0.969271 -0.0349723 -0.243511 0} {0 0 0 1}} {{0.0474491 0 0 0} {0 0.0474491 0 0} {0 0 0.0474491 0} {0 0 0 1}} {{1 0 0 -0.65} {0 1 0 -0.42} {0 0 1 0} {0 0 0 1}}}
	lappend viewplist [molinfo top]
	set topmol [molinfo top]
	# done with molecule 0
	foreach v $viewplist {
	  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)
	}
	foreach v $fixedlist {
	  molinfo $v set fixed 1
	}
	unset viewplist
	unset fixedlist
	#mol top $topmol
	unset topmol


}

proc rotateX {} {

	rotate x by 720 1
}

proc rotateY {} {

	rotate y by 720 1
}

proc rotateZ {} {

	rotate z by 720 1
}

proc zoomIn {} {

	scale by 1.2
}

proc zoomOut {} {

	scale by 0.833
}


proc Play {} {

	animate forward
}


proc Pause {} {

	animate pause
}


proc Rewind {} {

	animate reverse
}


proc animation {} {

	view4

	set nframes [molinfo top get numframes]

	for {set f 0} {$f <= $nframes } {incr f} {
		
		
		animate goto $f
		display update ui
	

		if {[molinfo top get frame] == 70} {rotate x by 180 2}

		if {[molinfo top get frame] == 160} {rotate y by 90 2}
	
		if {([molinfo top get frame] >= 71) && ([molinfo top get frame] <= 205) } {translate by 0 0 0.005}
		
		if {([molinfo top get frame] >= 206) && ([molinfo top get frame] <= 285) } {translate by 0 0 0.02}
	
		if {[molinfo top get frame] == 245} {rotate z by 180 4.5}
	
		if {[molinfo top get frame] == 285} {rotate y by -90 3}

		if {([molinfo top get frame] >= 286) && ([molinfo top get frame] <= 317) } {translate by 0 0 -0.05}
				
		if {([molinfo top get frame] >= 318) && ([molinfo top get frame] <= 344) } {translate by 0 0 -0.025}
				


	}

}


settings

