proc _bs {comp l key low high} {
  if {$low > $high} {return -1} else {
    set middle [expr {($low + $high)/2}]
    set test [$comp $key [lindex $l $middle]]
    if {$test == 0} {
      return $middle
    } elseif {$test < 0} {
      incr middle -1
      return [_bs $comp $l $key $low $middle] 
    } else {
      incr middle
      return [_bs $comp $l $key $middle $high]
    }
  }
}
 
proc blsearch {comp l key} {_bs $comp $l $key 0 [expr {[llength $l] - 1}]}
proc compint {x y} {expr {$x - $y}}