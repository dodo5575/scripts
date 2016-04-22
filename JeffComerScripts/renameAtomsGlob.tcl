source $env(HOME)/scripts/useful.tcl

set fileList [glob force*.pdb]
set name POT

foreach f $fileList {
    set base [trimExtension [trimExtension $f]]

    mol new $f
    set all [atomselect top all]
    $all set name $name
    $all writepdb $base.pdb
    
    $all delete
    mol delete top
}