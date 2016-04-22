#!/usr/bin/tclsh
# Make block averages.
# Use with: tclsh blockAverage.tcl window_size input_size
# jcomer2@uiuc.edu

if {$argc != 1} {
    puts "This script requires one argument: "
    puts "  input prefix"
    puts ""
    exit
}
proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}

# Input:
set prefix [trimExtension [lindex $argv 0]]
set inData ${prefix}.dat
set inLighting lighting.dat
# Output:
set outData ${prefix}_light.dat

# Get the lighting commands.
set in [open $inLighting r]
set lighting [read $in]
close $in

# Write the header.
set out [open $outData w]
# Get the data.
set in [open $inData r]
set found 0
while { [gets $in line] >= 0 } {
    set fields [split $line]
    case [lindex $fields 0] {
	Resolution {
	    set tok [concat $line]
	    set width [lindex $tok 1]
	    set height [lindex $tok 2]
	}
	Directional_Light {
	    if {!$found} {
		set line $lighting
		set found 1
	    } else { 
		set line ""
	    }
	}
    }
    puts $out $line
    if {$found > 0} {incr found}
    if {$found > 10} {
	break
    }
}
fcopy $in $out

close $in
close $out
set w [expr $width/4]
set h [expr $height/4]
puts "The tachyon file `${outData}' was written successfully."
puts "Test\ntachyon -fullshade -add_skylight 1.4 -rescale_lights 0.6 -res $w $h -aasamples 0 ${prefix}_light.dat -o ${prefix}_light.tga"
puts "Final\ntachyon -fullshade -add_skylight 1.4 -rescale_lights 0.6 -res $width $height -aasamples 4 ${prefix}_light.dat -o ${prefix}_light.tga"
