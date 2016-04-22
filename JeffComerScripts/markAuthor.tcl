# Add author information to scripts.
# Use with: tclsh markAuthor.tcl
# Author: Jeff Comer <jcomer2@illinois.edu>

if {$argc < 2} {
    puts "Usage: tclsh markAuthor.tcl outDir inputFile0 \[inputFile1\]..."
    exit
}

# Parameters:
set outDir [lindex $argv 0]
set inFileList [lrange $argv 1 end]

# Author: Jeff Comer <jcomer2@illinois.edu>
set authorLine "Author: Jeff Comer <jcomer2@illinois.edu>"

foreach fileName $inFileList {
    set ext [file extension $fileName]
    if {[string equal -nocase $ext ".tcl"]} {
	set comment "\#"
    } elseif {[string equal $ext ".C"]} {
	set comment "//"
    } elseif {[string equal -nocase $ext ".m"]} {
	set comment "%"
    } elseif {[string equal -nocase $ext ".pl"]} {
	set comment "\#"
    } else {
	puts "Unrecognized extension: $fileName"
	continue
    }

    set out [open $outDir/$fileName w]
    set in [open $fileName r]

    foreach line [split [read $in] \n] {
	
	if {[string match "*$comment*jcomer2@*" $line]} {
	    # Write the author line.
   	    puts $out "$comment $authorLine"
	} else {
	    # Just write any line that isn't an author line.
	    puts $out $line
	}
    }
    close $in
    close $out
}
exit



