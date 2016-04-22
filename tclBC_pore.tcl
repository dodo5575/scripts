#k = 1 kcal/(mol angstrom^2)
set atomName "P"
#set origamiForce [list 0 0 0]

#parsing the all atom pdb and get the list of "segname resid name" of selected atoms 
set origami_list {}

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
        lappend origami_list "[string trim $string2]"

    }

}
close $inStream

wrapmode input

proc calcforces {step unique k Zmin} {

    global origami_list radius 

    if {$step % 1200 == 0} {
	cleardrops
    }    

    set totalZForce 0

    while {[nextatom]} {

	set atomid [getid]

	#if {[lsearch $$origami_list $atomid] == -1} {
	#    
	#    dropatom
	#    continue
	#}

	foreach {x y z} [getcoord] {break}

	if {[lsearch $$origami_list $atomid] != -1 && $z < $Zmin && [expr $x*$x + $y*$y] < [expr $radius * $radius]} {

	    set mass [getmass]

	    set dZ [expr $Zmin - $z]
	    
	    set force [list 0 0 [expr $dZ * $k]]

	    set totalZForce [expr $totalZForce + [lindex $force 2]]

	    addforce $force

	    if {$dZ > 3} {print "Warning! Step:$step; Index:$atomid; Mass:$mass; Coordinate:($x,$y,$z) Displacement:$dZ; Force:$force (kcal/(mol angstrom)) "}
 
	} else {
	    dropatom

	}

    }

    if {$step % 1200 == 0} { 
        print "Step: $step ; Total z force: $totalZForce (kcal/(mol angstrom))"
    }
}
