set atomList {};
set atomID 1  ;
set sum 0 ;

set inStream [open conf.pdb r]   ;

foreach line [split [read $inStream] \n] {

  set string1 [string range $line 17 20]
  set string2 [string range $line 13 14]
  set res [string trim $string1]
  set atom [string trim $string2]
 
  print "atom index: $atomID $atom"
  if { [ string equal $atom {CA} ] } {
    #lappend atomList $atomID 
    lappend atomList 1
    print "atom index: $atomID $atom accepted"
  } else {
    lappend atomList 0
  }

  incr atomID
}
close $inStream

wrapmode cell

#foreach {x0 y0 z0} $bcCenter {break}

proc calcforces {step unique K} {

  global x0 y0 z0 droff bcCenter atomList sum

  foreach {x0 y0 z0 droff } $bcCenter { break }

  set zcore [expr $droff-5]

  while { [nextatom] } {

    set atomid [expr [getid]-1] ;

    if { [lindex $atomList $atomid] == 1 } {
        set atvec [getcoord] ;
        foreach {x y z} $atvec { break } ;

        set zdiff [expr $z-$z0]
        set forceZ 0
        set delZ 0
        if { $zdiff > $droff} {
          set delZ   [expr $zdiff-$droff]
          set forceZ [expr -$delZ*$K]
          #print "atomid: $atomid larger than $droff delZ: $delZ force: $forceZ"
	  addforce "0 0 $forceZ"
          set sum [expr  $sum + $forceZ]
        } elseif { $zdiff < -$droff} {
          set delZ   [expr $zdiff+$droff]
          set forceZ [expr -$delZ*$K]
          #print "atomid: $atomid less than -$droff delZ: $delZ force: $forceZ"
	  addforce "0 0 $forceZ"
          set sum [expr  $sum - $forceZ]
        } elseif { $zdiff < $zcore && $zdiff > -$zcore } {
	    #dropatom
	    #continue
        } else {
          set forceZ 0
        } 
      #print "step $step atom: $atomid zdiff: $zdiff  delz: $delZ forceZ: $forceZ constant: $K z: $z"
    } else {
	dropatom
    }
  }

  if { [expr $step % 400] == 0 } {
    set sum [expr $sum/400]
    print "$step : $sum"
    set sum 0
  }
}
