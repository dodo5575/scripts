######################################################################
# this is a tcl script for applying an oscillating electric field to
# selected atoms
# use with name tclforces
# written by alek@ks.uiuc.edu
# Author: Jeff Comer <jcomer2@illinois.edu>
######################################################################
# Inputs: eField0, forcesRecalcFreq, targetAtomPdb

print "Starting acField"

set pi [expr 4.*atan(1.)]
# define all protein atoms :
set targetMark    "1.00"

set prot_list {}
set chargeList {}
print "Reading the list of target atoms..."

# The beta column determines whether the force should be exerted and the
# occupancy column must contain the charge of the atom.
set inStream [open $targetAtomPdb r]
foreach line [split [read $inStream] \n] {
    set string1 [string range $line 0 3]
    set string2 [string range $line 6 10]
    set string3 [string range $line 17 20]	
    set string4 [string range $line 12 15]
    set string5 [string range $line 46 53]
    set string6 [string range $line 72 75]
    set string673 [string range $line 73 75]
    set string7 [string range $line 22 25]
    set string8 [string range $line 62 65]
    set string9 [string range $line 55 60]
   

    if { ([string equal $string1 {ATOM}] || \
	      [string equal $string1 {HETA}] ) && \
		 [string equal $targetMark $string8] } {
		
		set q [string trim $string9]
		if {![string equal $q "0.00"]} {
			lappend prot_list "[string trim $string6]\
			    [string trim $string7] [string trim $string4]"
			lappend chargeList $q
		}
    }  

}
close $inStream

# make list of indices


set protAtoms {}

foreach atomrecord $prot_list {
#    print $atomrecord
    foreach {segname resid atom} $atomrecord  { break }
    set atomindex [atomid $segname $resid $atom]
    lappend protAtoms $atomindex
    addatom  $atomindex
    print "$atomindex $segname $resid $atom"
}


set protAtoms   [concat $protAtoms]

# -- zero out initial forces --
set pushAtoms $protAtoms
set f {}
foreach index $pushAtoms {
	lappend f [list 0.0 0.0 0.0]
}

if { [llength $protAtoms] > 0 } {
    set push 1
} else {
    print "WARNING: no atoms have been selected"
    set push 0
}

print "Registered [llength $protAtoms] target atoms"


# initialize printing counter (independent on step counter)

set  pushCount [expr $forcesRecalcFreq -0]
set stepCount 0

###################################################################
# this procedure is executed at each time step
###################################################################

print "Starting calcforces..."
proc calcforces {} {

    global stepCount forcesRecalcFreq
    global pushCount stepCount
    global protAtoms pushAtoms
    global push 
    global f
    global eField0
    global chargeList
    global fieldPeriod
    global pi
      
     if {$push == 1} {


	##-------------  apply forces  --------------###	
	foreach i $pushAtoms force $f {
	    addforce $i $force
	    #	    print "atom $i force $force"
	}

	###------------ recalculate forces -----------###

	if { $pushCount == $forcesRecalcFreq } {		
    		##----- reconfiguting densities ------##
	
	    print "recalculating forces at $stepCount"	   

	    set eField1 [expr $eField0*sin(2.*$pi*$stepCount/$fieldPeriod)]
 		
	    set f {}
	    set pushAtoms $protAtoms
	    set nForce 0

	    foreach index $protAtoms charge $chargeList {	
		lappend f [list 0. 0. [expr $charge*$eField1]]
	    }
	    
	    print "step ${stepCount}: Recalculated [llength $f] forces "
	    set pushCount 0				
	}
	
	incr pushCount	 
    }
    incr stepCount
    return    
    
}



