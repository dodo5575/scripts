#!/usr/bin/tclsh

# Author: jcomer2@illinois.edu
set inFile doRemove.txt
set analysisScript analRemoveWater.tcl

set stride 1
set moiety nw
set outDir ../dcd

set struct(pos) ../extra_layer_pos
set struct(neg) ../extra_layer_neg
set struct(zer) ../extra_layer_zer
set struct(nng) ../extra_layer_nng
set struct(ngn) ../extra_layer_ngn
set struct(p2e) ../extra_layer_p2e
set struct(p2f) ../extra_layer_p2f

set struct(q-1) ../anneal_dopc_q-1
set struct(q0) ../anneal_dopc_q0
set struct(q1) ../anneal_dopc_q1
set struct(q2c) ../anneal_dopc_q2c
set struct(aq-1) ../anneal_dopc_q-1
set struct(aq0) ../anneal_dopc_q0
set struct(aq1) ../anneal_dopc_q1
set struct(aq2c) ../anneal_dopc_q2c

set struct(q2r) ../anneal_dopc_q2r
set struct(q2s) ../anneal_dopc_q2s

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr $ind+1] end]
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

set inp [open $inFile r]
while {[gets $inp line] >= 0} {
    if {[string length $line] <= 1} { continue }
    if {[string match "\#*" $line]} { continue }
    
    # Get the system name.
    set chunk [regexp -inline {dopc_[\-a-z0-9]*.dcd} $line]
    set sys [lindex [split $chunk "._"] 1]
    
    set name [trimPath [trimExtension $line]]
    if {![info exists struct($sys)]} {
	puts "ERROR: Unknown system `$sys' from $line!"
	exit
    }
    set structPrefix $struct($sys)

    puts "vmd -dispdev text -e $analysisScript -args $name $moiety $structPrefix $outDir $stride $line"    
}
close $inp
