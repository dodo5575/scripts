# Author: jcomer2@illinois.edu
set pdbGlob ../init/init_pot_at_pos*.pdb
set outPre "ref"

set fileList [glob $pdbGlob]

proc trimPath {name} {
    set ind [string last "/" $name]
    return [string range $name [expr {$ind+1}] end]
}
proc extractPath {name} {
    set ind [string last "/" $name]
    return [string range $name 0 [expr {$ind-1}]]
}

foreach f $fileList {
    # Load.
    mol load pdb $f
    set all [atomselect top all]
    set fileName [trimPath $f]
    set filePath [extractPath $f]

    # Set the beta column to the occupancy column.
    set occList [$all get occupancy]
    $all set beta $occList

    # Write the result.
    $all writepdb $filePath/${outPre}_${fileName}

    $all delete
    mol delete top
}
exit
