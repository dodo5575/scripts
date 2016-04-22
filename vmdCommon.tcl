# Author: Jeff Comer <jcomer2@illinois.edu>
puts "\n[string repeat "-" 10] vmdCommon.tcl [string repeat "-" 10]"


proc saveLights {{number 0}} {
  set saveName _saveLights
  set fileName ${saveName}${number}.dat

  if {[file exists $fileName]} {
    puts -nonewline "File `${fileName}' exists. Overwrite?"
    flush stdout
    gets stdin reply
    if {![string match -nocase "y*" $reply]} { return -1 }
  }
    set out [open $fileName w]
    for {set i 0} {$i<4} {incr i} {
      puts $out "[light $i pos] [light $i status]"
    }
    close $out
    puts "Saved the current lights  to `$fileName'."
    return 1

}

proc loadLights {{buffer 0}} {
  set saveName _saveLights
  set fileName ${saveName}${buffer}.dat

  set in [open $fileName r]
  set l 0
  while { [gets $in line] != -1 } {
    set data [split $line]
    light $l pos [lrange $data 0 2]
    light $l [lindex $data 3]
    puts "[light $l pos] [light $l status]"
    incr l
  }
  close $in
  puts "Let there be light! (as in `$fileName')."
}

proc saveView {{buffer 0}} {
  set saveName _saveView
  set fileName ${saveName}${buffer}.dat

  set bigOlView [molinfo top get {center_matrix rotate_matrix scale_matrix global_matrix}]

  if {[file exists $fileName]} {
    puts -nonewline "File `${fileName}' exists. Overwrite?"
    flush stdout
    gets stdin reply
    if {![string match -nocase "y*" $reply]} {
      return -1
    }
  }

  set out [open $fileName w]
  puts $out $bigOlView
  close $out

  puts "Saved the current view to `$fileName'."

  return $bigOlView
}

proc loadView {{buffer 0}} {
  set saveName _saveView
  set fileName ${saveName}${buffer}.dat

  set in [open $fileName r]
  if {[gets $in line] < 0} {
    puts "Could not read `$fileName'."
    return -1
  }

  set bigOlView [concat $line]
  close $in

  puts "Loaded the view from `$fileName'."

  set molList [molinfo list]
  foreach m $molList {
    molinfo $m set {center_matrix rotate_matrix scale_matrix global_matrix} $bigOlView
  }
  return $bigOlView
}

proc setColors {} {
  color change rgb 0 0.2199999988079071 0.20999999344348907 0.8999999761581421
  color change rgb 2 0.3499999940395355 0.3499999940395355 0.3499999940395355
  color change rgb 3 0.36000001430511475 0.6299999952316284 0.25
  color change rgb 4 0.8799999952316284 0.8799999952316284 0.30000001192092896
  color change rgb 5 0.5 0.5 0.20000000298023224
  color change rgb 6 0.7300000190734863 0.7300000190734863 0.7300000190734863
  color change rgb 7 0.0 1.0 0.0
  color change rgb 9 0.8999999761581421 0.5199999809265137 0.5199999809265137
  color change rgb 11 0.6499999761581421 0.0 0.6499999761581421
  color change rgb 12 0.5 0.8999999761581421 0.4000000059604645
  color change rgb 13 0.8999999761581421 0.4000000059604645 0.699999988079071
  color change rgb 14 0.5 0.30000001192092896 0.0
  color change rgb 15 0.4699999988079071 0.4699999988079071 0.8999999761581421
  color change rgb 17 0.8999999761581421 0.9399999976158142 0.800000011920929
  color change rgb 18 0.550000011920929 0.8999999761581421 0.019999999552965164
  color change rgb 19 0.0 0.8999999761581421 0.03999999910593033
  color change rgb 20 0.0 0.8999999761581421 0.5
  color change rgb 21 0.0 0.8799999952316284 1.0
  color change rgb 22 0.0 0.6299999952316284 0.8899999856948853
  color change rgb 23 0.019999999552965164 0.3799999952316284 0.6700000166893005
  color change rgb 24 0.03999999910593033 0.10000000149011612 0.949999988079071
  color change rgb 25 0.27000001072883606 0.0 0.9800000190734863
  color change rgb 26 0.44999998807907104 0.0 0.8999999761581421
  color change rgb 27 0.8999999761581421 0.6200000047683716 0.8999999761581421
  color change rgb 28 0.6200000047683716 0.3100000023841858 1.0
  color change rgb 29 0.9800000190734863 0.0 0.23000000417232513
  color change rgb 30 0.8999999761581421 0.15000000596046448 0.15000000596046448
  color change rgb 31 0.8899999856948853 0.3499999940395355 0.0
  color change rgb 32 0.9599999785423279 0.7200000286102295 0.1599999964237213
}


proc setColors1 {} {
  color change rgb 0    0.00  0.00  1.00
  color change rgb 2    0.35  0.35  0.35
  color change rgb 3    1.00  0.50  0.00
  color change rgb 4    1.00  1.00  0.00
  color change rgb 5    0.50  0.50  0.20
  color change rgb 6    0.60  0.60  0.60
  color change rgb 7    0.00  1.00  0.00
  color change rgb 9    1.00  0.60  0.60
  color change rgb 11   0.65  0.00  0.65
  color change rgb 12   0.50  0.90  0.40
  color change rgb 13   0.90  0.40  0.70
  color change rgb 14   0.50  0.30  0.00
  color change rgb 15   0.50  0.50  0.75
  color change rgb 17   0.88  0.97  0.02
  color change rgb 18   0.55  0.90  0.02
  color change rgb 19   0.00  0.90  0.04
  color change rgb 20   0.00  0.90  0.50
  color change rgb 21   0.00  0.88  1.00
  color change rgb 22   0.00  0.76  1.00
  color change rgb 23   0.02  0.38  0.67
  color change rgb 24   0.01  0.04  0.93
  color change rgb 25   0.27  0.00  0.98
  color change rgb 26   0.45  0.00  0.90
  color change rgb 27   0.90  0.00  0.90
  color change rgb 28   1.00  0.00  0.66
  color change rgb 29   0.98  0.00  0.23
  color change rgb 30   0.81  0.00  0.00
  color change rgb 31   0.89  0.35  0.00
  color change rgb 32   0.96  0.72  0.00
}

proc setColors2 {} {
  color change rgb 0 0.2199999988079071 0.20999999344348907 0.8999999761581421
  color change rgb 2 0.3499999940395355 0.3499999940395355 0.3499999940395355
  color change rgb 3 0.36000001430511475 0.6299999952316284 0.25
  color change rgb 4 0.8799999952316284 0.8799999952316284 0.30000001192092896
  color change rgb 5 0.5 0.5 0.20000000298023224
  color change rgb 6 0.7300000190734863 0.7300000190734863 0.7300000190734863
  color change rgb 7 0.0 1.0 0.0
  color change rgb 9 0.8999999761581421 0.5199999809265137 0.5199999809265137
  color change rgb 11 0.6499999761581421 0.0 0.6499999761581421
  color change rgb 12 0.5 0.8999999761581421 0.4000000059604645
  color change rgb 13 0.8999999761581421 0.4000000059604645 0.699999988079071
  color change rgb 14 0.5 0.30000001192092896 0.0
  color change rgb 15 0.4699999988079071 0.4699999988079071 0.8999999761581421
  color change rgb 17 0.8799999952316284 0.9700000286102295 0.019999999552965164
  color change rgb 18 0.550000011920929 0.8999999761581421 0.019999999552965164
  color change rgb 19 0.16 0.66 0.27
  color change rgb 20 0.0 0.8999999761581421 0.5
  color change rgb 21 0.0 0.8799999952316284 1.0
  color change rgb 22 0.0 0.6299999952316284 0.8899999856948853
  color change rgb 23 0.019999999552965164 0.3799999952316284 0.6700000166893005
  color change rgb 24 0.03999999910593033 0.10000000149011612 0.949999988079071
  color change rgb 25 0.27000001072883606 0.0 0.9800000190734863
  color change rgb 26 0.44999998807907104 0.0 0.8999999761581421
  color change rgb 27 0.8999999761581421 0.6200000047683716 0.8999999761581421
  color change rgb 28 0.6200000047683716 0.3100000023841858 1.0
  color change rgb 29 0.9800000190734863 0.0 0.23000000417232513
  color change rgb 30 0.8999999761581421 0.15000000596046448 0.15000000596046448
  color change rgb 31 0.8899999856948853 0.3499999940395355 0.0
  color change rgb 32 0.9599999785423279 0.7200000286102295 0.1599999964237213
}

proc representNucleic {} {
  set n [molinfo top get numreps]
  for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

  display depthcue on
  light 0 on
  mol selection {nucleic}
  mol color Name
  mol representation Licorice
  mol material Opaque
  mol addrep top
}

proc representSin {} {
  set n [molinfo top get numreps]
  for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

  display projection orthographic
  display resetview
  rotate y by 90
  color Display {Background} white
  display depthcue off
  light 0 on
  light 1 off
  mol representation MSMS
  mol color ColorID 0
  mol selection {resname SIN}
  mol material Transparent
  mol addrep top
  mol drawframes top 0 {0}
}

proc representPore {} {
  representSin

  mol representation VDW
  mol color Name
  mol selection {ions}
  mol material Opaque
  mol addrep top

  ruler grid
}


proc representVolSlice {{center 0.0} {scale 6.0}} {
  set n [molinfo top get numreps]
  for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

  display depthcue off
  light 0 off
  color Display {Background} gray
  display backgroundgradient off
  color Display {FPS} white
  color Axes {Labels} white
  color Labels {Bonds} white
  color Labels {Angles} blue
  color scale method BWR

  mol representation VolumeSlice 0.500000 0.000000 2.000000 0.000000
  mol color Volume 0
  mol selection {all}
  mol material Opaque
  mol addrep top

  mol scaleminmax top 0 [expr {$center-0.5*$scale}] [expr {$center+0.5*$scale}]
}

proc setRadii {} {
  set pot [atomselect top "name POT"]
  $pot set radius 1.76375
  $pot delete

  set cla [atomselect top "name CLA"]
  $cla set radius 2.27
  $cla delete

  set si [atomselect top "resname SIO2 SIO and name \"SI.*\""]
  $si set radius 2.1475
  $si delete

  set osi [atomselect top "resname SIO2 SIO and name \"O.*\""]
  $osi set radius 1.75
  $osi delete
}

proc representDnaSpheres {{mat Opaque} {radius 1.0}} {
  set backbone "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"

  #set n [molinfo top get numreps]
  #for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

  set nuc [atomselect top nucleic]
  set segList [lsort -unique [$nuc get segname]]
  $nuc delete

  set darkList {30 0}
  set lightList {9 15}

  set j 0
  foreach seg $segList {

    set dark [lindex $darkList [expr {$j % [llength $darkList]}]]
    set light [lindex $lightList [expr {$j % [llength $lightList]}]]

    mol representation VDW $radius 15.000000
    mol color ColorID $dark
    mol selection "segname $seg and name \"C.*\" \"O.*\""
    mol material $mat
    mol addrep top

    mol representation VDW $radius 15.000000
    mol color ColorID $light
    mol selection "segname $seg and not name \"C.*\" \"O.*\""
    mol material $mat
    mol addrep top

    incr j
  }

}

proc representDnaBackbone {{mat Opaque} {radius 1.0}} {
  set backbone "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"

  #set n [molinfo top get numreps]
  #for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

  set nuc [atomselect top nucleic]
  set segList [lsort -unique [$nuc get segname]]
  $nuc delete

  set darkList {30 0}
  set lightList {9 15}

  set j 0
  foreach seg $segList {

    set dark [lindex $darkList [expr {$j % [llength $darkList]}]]
    set light [lindex $lightList [expr {$j % [llength $lightList]}]]

    mol representation VDW $radius 15.000000
    mol color ColorID $dark
    mol selection "segname $seg and not name $backbone"
    mol material $mat
    mol addrep top

    mol representation VDW $radius 15.000000
    mol color ColorID $light
    mol selection "segname $seg and name $backbone"
    mol material $mat
    mol addrep top

    incr j
  }

}

proc representDnaCartoon {{mat Opaque}} {
  set backbone "C1' H1' C2' H2' H2'' C3' O3' H3' C4' O4' H4' C5' O5' H5' H5'' O1P O2P P"

  set n [molinfo top get numreps]
  for {set i 0} {$i < $n} {incr i} { mol delrep 0 top }

  set nuc [atomselect top nucleic]
  set segList [lsort -unique [$nuc get segname]]
  $nuc delete

  set darkList {30 0}
  set lightList {9 15}

  set j 0
  set count 0
  foreach seg $segList {

    set dark [lindex $darkList [expr {$j % [llength $darkList]}]]
    set light [lindex $lightList [expr {$j % [llength $lightList]}]]

    mol representation NewCartoon 0.900000 20.000000 1.280000 0
    mol color ColorID $dark
    mol selection "segname $seg"
    mol material $mat
    mol addrep top
    incr count

    mol representation VDW 0.8 15.000000
    mol color ColorID $dark
    mol selection "segname $seg and name C3'"
    mol material $mat
    mol addrep top
    incr count

    incr j
  }

  for {set i 0} {$i < $count} {incr i} {
    mol smoothrep top $i 0
  }
}


proc representIons {{mat Opaque}} {
  setRadii

  mol representation VDW 1.000000 15.000000
  mol color ColorID 17
  mol selection "name POT SOD"
  mol material $mat
  mol addrep top

  mol representation VDW 1.000000 15.000000
  mol color ColorID 18
  mol selection "name CLA"
  mol material $mat
  mol addrep top
}

proc silent {} {
}
puts "\n[string repeat "-" 3] Finished with vmdCommon.tcl [string repeat "-" 3]"


