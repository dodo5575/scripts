

proc assignSeg {inPut outPut} {

	set in [open $inPut.pdb r]		
	set out [open $outPut.pdb w]

	while {[gets $in line] >= 0} {

		if {[lindex $line 0] == "ATOM"}	{

			set atomNum [lindex $line 1]
	
			if {$atomNum <= 660} {
				set line "$line\tADNA"
				
			}
			if {$atomNum >= 661 && $atomNum <= 1260} {
				set line "$line\tBDNA"
				
			}
			
			
			
		}
		
		set resName [lindex $line 3]
		
		if {$resName == "NA"} {
			set line "$line\tSOD"
		}
		if {$resName == "CL"} {
			set line "$line\tCLA"
		}
		if {$resName == "SOL"} {
			set line "$line\tWAT"
		}
		

		puts $out $line 

	}
	close $in
	close $out

}




