#!/usr/bin/tclsh
# Make block averages.
# Use with: tclsh blockAverage.tcl window_size input_size
# jcomer2@uiuc.edu

if {$argc != 1} {
	puts "This script requires one argument: "
	puts "  input file"
	puts ""
	exit
}

proc trimExtension {name} {
    set ind [string last "." $name]
    return [string range $name 0 [expr $ind-1]]
}


# Input:
set prefix [trimExtension [lindex $argv 0]]
set inData ${prefix}.pov
set inLighting lighting.pov
# Output:
set outData ${prefix}_light.pov

# Get the lighting commands.
set in [open $inLighting r]
set lighting [read $in]
close $in

# Write the header.
set out [open $outData w]
# Get the data.
set in [open $inData r]
set onLight -1
foreach line [split [read $in] \n] {
	if {[string match "light_source*" $line]} {
		if {$onLight < 0} {
			puts $out $lighting
			set onLight 1
		} else {
			set onLight 1
		}
		puts "Begin light_source block"
	} elseif {$onLight < 1} {	
		if {[string match "// try povray*" $line]} {
			set param [split $line]
			set width [string range [lindex $param 3] 2 end]
			set height [string range [lindex $param 4] 2 end]
			puts "Dimensions: $width $height"
			puts $out "// try povray +W${width} +H${height} -I${outData} -O${outData}.tga +P +X +A +FT"
		} else {
			puts $out $line
		}
	} elseif {[string match "*\}*" $line]} {
		set onLight 0
		puts "End light_source block"
	}
}
close $in


puts "The povray file `${outData}' was written successfully."
puts "try `povray +W${width} +H${height} -I${outData} -O${outData}.tga +P +X +A +FT'"
close $out


