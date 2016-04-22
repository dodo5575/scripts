
puts "draw_exb_atom input"

proc draw_exb_atom {input} {

    set indices ""
    
    set inStream [open $input r]
    foreach line [split [read $inStream] \n] {
    
        regexp {bond\s*([0-9]*)\s*([0-9]*)} $line matched ind1 ind2
        lappend indices $ind1
        lappend indices $ind2
    
    }
    close $inStream
    
    set indices [lsort -integer -unique $indices]
    
    puts $indices 
    
    set NumOfRep [lindex [mol list top] 12]
    
    mol addrep top
    mol modcolor    $NumOfRep top ColorID 16
    mol modstyle    $NumOfRep top VDW 1.000000 12.000000
    mol modselect   $NumOfRep top index $indices
    mol modmaterial $NumOfRep top BrushedMetal

}

