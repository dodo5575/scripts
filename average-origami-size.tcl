
proc averageOrigamiSize {PSF dimension outPut dcdlist} {

        set outFile [open $outPut a+]

	mol new $PSF
        foreach dcdFile $dcdlist {
	    mol addfile $dcdFile waitfor all
        }


	set DNA [atomselect top nucleic]
	set nFrames [molinfo top get numframes]
	
	set length 0
	
	for {set f 0} {$f <= [expr $nFrames - 1]} {incr f} {
	$DNA frame $f
	set l [measure minmax $DNA]
	
	switch -- $dimension {
		"x" {
			set dim 0
        	}
        	"y" {
        		set dim 1
        	}
        	"z" {
        		set dim 2
        	}
        	default {
        	}
	}

	set L1 [lindex $l 0 $dim] 
	set L2 [lindex $l 1 $dim]
	
	set length [expr {abs($L1 - $L2) + $length} ]
	}
	
	puts $outFile [expr {double($length) / $nFrames}]

        mol delete all
        close $outFile
}

if {$argc < 1} {

	puts "vmd -dispdev text -e $argv0 -args PSF dimension outPut dcdlist" 
}

set PSF [lindex $argv 0]
#set DCD [lindex $argv 1]
set dimension [lindex $argv 1]
set outPut [lindex $argv 2]
set dcdList [lrange $argv 3 end]



averageOrigamiSize $PSF $dimension $outPut $dcdList

exit
