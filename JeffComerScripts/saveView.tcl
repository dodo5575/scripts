# Author: Jeff Comer <jcomer2@illinois.edu>

proc saveView {{buffer 0}} {
    set saveName _saveView
    set fileName ${saveName}${buffer}.dat

    set bigOlView [molinfo top get {center_matrix rotate_matrix scale_matrix global_matrix}]

    if {[file exists $fileName]} {
	puts -nonewline "File `${fileName}' exists. Overwrite?"
	flush stdout
	gets stdin reply
	if {![string match -nocase "y*" $reply]} {
	    return -1
	}
    }

    set out [open $fileName w]
    puts $out $bigOlView
    close $out

    puts "Saved the current view to `$fileName'."

    return $bigOlView
}

proc loadView {{buffer 0}} {
    set saveName _saveView
    set fileName ${saveName}${buffer}.dat
    
    set in [open $fileName r]
    if {[gets $in line] < 0} {
	puts "Could not read `$fileName'."
	return -1
    }

    set bigOlView [concat $line]
    close $in

    puts "Loaded the view from `$fileName'."

    set molList [molinfo list]
    foreach m $molList {
	molinfo $m set {center_matrix rotate_matrix scale_matrix global_matrix} $bigOlView
    }
    return $bigOlView
}
