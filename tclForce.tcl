
#k = 1 kcal/(mol angstrom^2)
set atomName "P"
set origamiForce [list 0 0 0]
set sampleDNAForce [list 0 0 0]

#parsing the all atom pdb and get the list of "segname resid name" of selected atoms 
set origami_list {}
set sampleDNA_list {}

set inStream [open $allatompdb r]
foreach line [split [read $inStream] \n] {

    set string1 [string range $line 0 3]	
    #ATOM

    set string2 [string range $line 6 10]	
    #index 

    set string3 [string trim [string range $line 12 15]] 	
    #name

    set string4 [string range $line 17 19] 	
    #resname

    set string5 [string range $line 22 25]      
    #resid

    set string6 [string range $line 46 53] 	
    #z-coordinate

    set string7 [string range $line 72 75] 	
    #segname

#print "name $string3"
#print "resid $string5"
#print "segname $string7"

    if {[string equal {ATOM} $string1] && \
	(([string match SC* $string7] || \
	[string match P* $string7] ) && \
	[string equal $atomName $string3])} {
	#print "true"
	lappend origami_list "[string trim $string7] [string trim $string5] $atomName"

    }

    if {[string equal {ATOM} $string1] && \
	([string match DS* $string7] && \
	[string equal $atomName $string3])} {

	lappend sampleDNA_list "[string trim $string7]\
				[string trim $string5] $atomName"

    }
}
close $inStream

#print "origami_list $origami_list"

#make lisf of indices
set origami {}
set sampleDNA {}

foreach atomrecord $origami_list {
    foreach {segname resid atom} $atomrecord {break}
    set atomindex [atomid $segname $resid $atom]
    lappend origami $atomindex
    addatom $atomindex
}

foreach atomrecord $sampleDNA_list {
    foreach {segname resid atom} $atomrecord {break}
    set atomindex [atomid $segname $resid $atom]
    lappend sampleDNA $atomindex
    addatom $atomindex
}


set count 0

proc calcforces { } {

    global kZ kXY
    global origami sampleDNA
    global count checkCount
    global coords_origami forces_origami
    global coords_sampleDNA forces_sampleDNA
    global origamiForce sampleDNAForce 
    global origamiX origamiY origamiZ sampleDNAX sampleDNAY sampleDNAZ 


    if { 1 } {
        # Apply forces
        foreach atom $origami {
            addforce $atom $origamiForce
        }

        foreach atom $sampleDNA {
            addforce $atom $sampleDNAForce
        }



	if {[expr $count % $checkCount] == 0} {
		
	    foreach atom $origami {
		addatom $atom
	    }
	}
	
	if {[expr $count % $checkCount] == 1} {	

		array unset coords_origami
		loadcoords coords_origami

		set SIZE [array size coords_origami]
		#print "$count\tcoords_origami array $SIZE"

		set origamiX 0
		foreach index $origami {
		    set origamiX [expr $origamiX + [lindex $coords_origami($index) 0]]
		}
		
		set origamiX [expr $origamiX / $SIZE]

		set origamiY 0
		foreach index $origami {
		    set origamiY [expr $origamiY + [lindex $coords_origami($index) 1]]
		}
		
		set origamiY [expr $origamiY / $SIZE]

		set origamiZ 0
		foreach index $origami {
		    set origamiZ [expr $origamiZ + [lindex $coords_origami($index) 2]]
		}
		
		set origamiZ [expr $origamiZ / $SIZE]


		#print "origamiZ $origamiZ"

		if {$count > $checkCount} {


			set dx [expr $sampleDNAX - $origamiX]
			set dy [expr $sampleDNAY - $origamiY]
			set dz [expr $sampleDNAZ - $origamiZ]
			

			set ForceX [expr $kXY * $dx / $SIZE]
			set ForceY [expr $kXY * $dy / $SIZE]
			set ForceZ [expr $kZ * $dz / $SIZE]

			set origamiForce [list $ForceX $ForceY $ForceZ]
			print "$count\tcenter of origami = $origamiX $origamiY $origamiZ"
			print "$count\tcenter of sampleDNAZ = $sampleDNAX $sampleDNAY $sampleDNAZ"
			print "$count\tcalculated force of origami : $origamiForce (kcal/(mol ang))"
		
		}

	}

	if {($count > $checkCount) && [expr $count % $checkCount] == 3} {

		array unset forces_origami
		loadforces forces_origami

		set SIZE [array size forces_origami]
		set NAMES [array names forces_origami]

		#print "$count\tforces_origami array $SIZE"
		#print "forces_origami array index $NAMES"
		#foreach index $origami {
		#
		#    print "measured origami forces $index : $forces_origami($index)"		
	
		#}

	}

	if {[expr $count % $checkCount] == 4} {
	
	    foreach atom $sampleDNA {
		addatom $atom
	    }
	}

	if {[expr $count % $checkCount] == 5} {
		array unset coords_sampleDNA
		loadcoords coords_sampleDNA

		set SIZE [array size coords_sampleDNA]
		#print "$count\tcoords_sampleDNA array $SIZE"

		set sampleDNAX 0
		foreach index $sampleDNA {
		    set sampleDNAX [expr $sampleDNAX + [lindex $coords_sampleDNA($index) 0]]
		}
		set sampleDNAX [expr $sampleDNAX / $SIZE]

		set sampleDNAY 0
		foreach index $sampleDNA {
		    set sampleDNAY [expr $sampleDNAY + [lindex $coords_sampleDNA($index) 1]]
		}
		set sampleDNAY [expr $sampleDNAY / $SIZE]

		set sampleDNAZ 0
		foreach index $sampleDNA {
		    set sampleDNAZ [expr $sampleDNAZ + [lindex $coords_sampleDNA($index) 2]]
		}
		set sampleDNAZ [expr $sampleDNAZ / $SIZE]

		if {$count > $checkCount} {
		

			set dx [expr $origamiX - $sampleDNAX]
			set dy [expr $origamiY - $sampleDNAY]
			set dz [expr $origamiZ - $sampleDNAZ]

			set ForceX [expr $kXY * $dx / $SIZE]
			set ForceY [expr $kXY * $dy / $SIZE]
			set ForceZ [expr $kZ * $dz / $SIZE]

			set sampleDNAForce [list $ForceX $ForceY $ForceZ]

			print "$count\tcenter of origami = $origamiX $origamiY $origamiZ"
			print "$count\tcenter of sampleDNA = $sampleDNAX $sampleDNAY $sampleDNAZ"
			print "$count\tcalculated force of sampleDNA : $sampleDNAForce (kcal/(mol ang))"


		}

	}

	if {($count > $checkCount) && [expr $count % $checkCount] == 7} {

		array unset forces_sampleDNA
		loadforces forces_sampleDNA

		set SIZE [array size forces_sampleDNA]
		set NAMES [array names forces_sampleDNA]

		#print "$count\tforces_sampleDNA array $SIZE"
		#print "forces_sampleDNA array index $NAMES"
		#foreach index $sampleDNA {
		#
		#    print "measured sampleDNA forces $index : $forces_sampleDNA($index)"		
	
		#}

	}


	if {[expr $count % $checkCount == 3] || [expr $count % $checkCount == 7]} {clearconfig}


    }
    incr count
}


