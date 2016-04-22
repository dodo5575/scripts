# Make extrabonds files for H-bond from pdb
# Usage: tclsh Hbond.tcl inPDB out 
# Chen-Yu Li	cli56@illinois.edu
# 2014/1/14


#pre-defined function
proc vecSub {a b} {
    set c {}
    foreach ai $a bi $b {
        lappend c [expr $ai-$bi]
    }
    return $c
}

proc vecLength {a} {
    set sum 0
    foreach ai $a {
        set sum [expr $sum + $ai*$ai]
    }
    return [expr sqrt($sum)]
}

# set variables and lists
set k 50

set inPDB [lindex $argv 0]
set out [lindex $argv 1]

set inStream [open $inPDB r]
set outFile [open $out w]

set G_O6 {}
set G_H1 {}
set G_H21 {}

set C_H42 {}
set C_N3 {}
set C_O2 {}

set A_H62 {}
set A_N1 {}

set T_O4 {}
set T_H3 {}


#parsing the all atom pdb and get the list of "segname resid name" of selected atoms 

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

    set string6 [string range $line 30 37]
    #x-coordinate

    set string7 [string range $line 38 45]
    #y-coordinate

    set string8 [string range $line 46 53]
    #z-coordinate

    set string9 [string range $line 72 75]
    #segname

#print "name $string3"
#print "resid $string5"
#print "segname $string7"

    if {[string equal {ATOM} $string1] && ([string match GUA $string4] && [string equal {O6} $string3])} {
        #print "true"
        lappend G_O6 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }
    if {[string equal {ATOM} $string1] && ([string match GUA $string4] && [string equal {H1} $string3])} {
        #print "true"
        lappend G_H1 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }
    if {[string equal {ATOM} $string1] && ([string match GUA $string4] && [string equal {H21} $string3])} {
        #print "true"
        lappend G_H21 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }
    if {[string equal {ATOM} $string1] && ([string match CYT $string4] && [string equal {H42} $string3])} {
        #print "true"
        lappend C_H42 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }
    if {[string equal {ATOM} $string1] && ([string match CYT $string4] && [string equal {N3} $string3])} {
        #print "true"
        lappend C_N3 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }
    if {[string equal {ATOM} $string1] && ([string match CYT $string4] && [string equal {O2} $string3])} {
        #print "true"
        lappend C_O2 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }
    if {[string equal {ATOM} $string1] && ([string match ADE $string4] && [string equal {H62} $string3])} {
        #print "true"
        lappend A_H62 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }
    if {[string equal {ATOM} $string1] && ([string match ADE $string4] && [string equal {N1} $string3])} {
        #print "true"
        lappend A_N1 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }
    if {[string equal {ATOM} $string1] && ([string match THY $string4] && [string equal {O4} $string3])} {
        #print "true"
        lappend T_O4 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }
    if {[string equal {ATOM} $string1] && ([string match THY $string4] && [string equal {H3} $string3])} {
        #print "true"
        lappend T_H3 "[string trim $string2] [string trim $string3] [string trim $string4] [string trim $string5] [string trim $string6] [string trim $string7] [string trim $string8] [string trim $string9]"

    }

}

close $inStream



#index name resname resid X Y Z segname
foreach GO6 $G_O6 {

    foreach CH42 $C_H42 {
    
	set vec1 [lrange $GO6 4 6]
	set vec2 [lrange $CH42 4 6] 

	set index1 [lindex $GO6 0]
	set index2 [lindex $CH42 0]

	set distance [vecLength [vecSub $vec1 $vec2]] 

	if {$distance <= 2.9} {

	    puts $outFile [format "bond%10d%10d%10.3g%10.3g" [expr $index1 - 1] [expr $index2 -1] $k 1.975]

	}

    }

}

#index name resname resid X Y Z segname
foreach GH1 $G_H1 {

    foreach CN3 $C_N3 {
    
	set vec1 [lrange $GH1 4 6]
	set vec2 [lrange $CN3 4 6] 

	set index1 [lindex $GH1 0]
	set index2 [lindex $CN3 0]

	set distance [vecLength [vecSub $vec1 $vec2]] 

	if {$distance <= 2.9} {

	    puts $outFile [format "bond%10d%10d%10.3g%10.3g" [expr $index1 - 1] [expr $index2 -1] $k 1.965]

	}

    }

}

#index name resname resid X Y Z segname
foreach GH21 $G_H21 {

    foreach CO2 $C_O2 {
    
	set vec1 [lrange $GH21 4 6]
	set vec2 [lrange $CO2 4 6] 

	set index1 [lindex $GH21 0]
	set index2 [lindex $CO2 0]

	set distance [vecLength [vecSub $vec1 $vec2]] 

	if {$distance <= 2.9} {

	    puts $outFile [format "bond%10d%10d%10.3g%10.3g" [expr $index1 - 1] [expr $index2 -1] $k 1.885]

	}

    }

}

#index name resname resid X Y Z segname
foreach AH62 $A_H62 {

    foreach TO4 $T_O4 {
    
	set vec1 [lrange $AH62 4 6]
	set vec2 [lrange $TO4 4 6] 

	set index1 [lindex $AH62 0]
	set index2 [lindex $TO4 0]

	set distance [vecLength [vecSub $vec1 $vec2]] 

	if {$distance <= 2.9} {

	    puts $outFile [format "bond%10d%10d%10.3g%10.3g" [expr $index1 - 1] [expr $index2 -1] $k 1.885]

	}

    }

}

#index name resname resid X Y Z segname
foreach AN1 $A_N1 {

    foreach TH3 $T_H3 {
    
	set vec1 [lrange $AN1 4 6]
	set vec2 [lrange $TH3 4 6] 

	set index1 [lindex $AN1 0]
	set index2 [lindex $TH3 0]

	set distance [vecLength [vecSub $vec1 $vec2]] 

	if {$distance <= 2.9} {

	    puts $outFile [format "bond%10d%10d%10.3g%10.3g" [expr $index1 - 1] [expr $index2 -1] $k 1.915]

	}

    }

}

