# This script generates the data mining index file used
# the Matlab script.
# It alerts you if any files fail to exist.
# Author: Jeff Comer <jcomer2@illinois.edu>

# Output:
set outFile index_curv.txt

# Defines:
set shift .shifted

set out [open $outFile w]
puts $out "# Data mining index file"
puts $out "# Format: "
puts $out "# name long_name voltage COM_position current position stretch curvature"

set num 0
set dir [list ../center_of_mass ../carrier ../pos ../stretch ../curvature]
set prefix [list com_ curr_ "" str_ curv_]
set extension [list .dat .dat .txt .dat .dat]

# Output a line for the index file.
# %n will be replaced with the fileName.
# %p uses the usual prefix.
# %x uses the usual extension.
# %u is equivalent to %p%n%x.
# %s is ${shift}
proc add {name longName voltage files} {
    global out dir num prefix extension shift
    
    # Write identifiers.
    puts -nonewline $out "$name $voltage"
    puts $name

    # Write the files.
    foreach f $files d $dir p $prefix x $extension {
	# Parse the file text, making substitutions.
	set n [string length $f]
	set s ""
	for {set j 0} {$j < $n} {incr j} {
	    # Get the current character.
	    set c [string index $f $j]

	    if {[string equal $c "%"]} {
		# Perform the substitutions.
		set cNext [string index $f [expr $j+1]]
		switch $cNext {
		    n {set s ${s}${longName}}
		    p {set s ${s}${p}}
		    x {set s ${s}${x}}
		    u {set s ${s}${p}${longName}${x}}
		    s {set s ${s}${shift}}
		}
		# Skip the next character.
		incr j
		
	    } else {
		# Just add this character.
		set s ${s}${c}
	    }
	}

	# Check that the file exists.
	set fileName $d/$s
	if {![file exists $fileName]} {
	    puts stderr "WARNING: $fileName does not exist."
	}
	
	# Write the file.
	puts -nonewline $out " ${fileName}"
    }
    puts $out ""
    incr num
}

# Begin the file indices.
# COM_position current position stretch curvature
add loop_2V loop_first_2V0-3 2.0 {%u %u %u %u %u}
add loop_1.5V loop_first_1.5V0-9 1.5 {%u %ploop_1.5V0-9%x %u %u %u}
add loop_6Vo hairpin_first_6V1-2 6.0 {%u %u head_6V1-2%x %u %u}
add loop_4Vo hairpin_first_4V0-2 4.0 {%u %u head_4V0-2%x %u %u}
add loop_3Vo hairpin_first_3V0-2 3.0 {%u %u head_3V0-2%x %u %u}
add loop_2.5Vo hairpin_first_2.5V0-6 2.5 {%u %u head_2.5V0-6%x %u %u}
add loop_1Vo hairpin_first_1V0-9 1.0 {%u %u %u %u %u}
add coil_2V coil_first_2V0-7 2.0 {%u %u %u %u %u}
add coil_1.5V coil_first_1V0-1.5V4 1.5 {%u %u %u %u %u}
add coil_cont_1.5V coil_cont_1.5V0-6 1.5 {%u%s %pcoil_first_1.5V0-6%x%s %u%s %u%s %u%s}
add coil_3V coil_first_3V0-2 3.0 {%u %u coil_first_3V%x %u %u}
add coil_4V coil_first_4V-5 4.0 {%u %pcoil_first_4V-6%x coil_first_4Va%x %u %u}
add coil_4Vo coil_first_E4V1-9 4.0 {%u %pE4V_1-9.txt coil_first_E4V%x %u %u}
add coil_6V coil_first_6V1-5 6.0 {%u%s %u%s coil_first_E6V_2ns%x%s %u%s %u%s}
if {0} {
add double_4V double_4V0-2 4.0 {%u %u %u %u %u}
add double_4Vo DS4_4V1-5 4.0 {%u %pDS4_4V.txt pos_DS4_4V1-5%x %u %u}
add double_5Vo double_5V1-4 5.0 {%u %pDS4_5V.txt pos_DS4_5V1-4%x %u %u}
add double_6Vo double_6V1-4 6.0 {%u %u pos_DS4_6V1-4%x %u %u}
}
close $out



