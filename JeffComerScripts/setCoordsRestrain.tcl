foreach sys {neg pos neu } {
    set selText "resname SIO2"
    set beta 20.0
    # Input:
    set refName pore40+dopc_${sys}1
    set fitName extra_layer_${sys}
    # Output:
    set outName layer_restrain_${sys}

    set refMol [mol load psf $refName.psf pdb $refName.pdb]
    set fitMol [mol load psf $fitName.psf pdb $fitName.pdb]

    set refSel [atomselect $refMol $selText]
    set fitSel [atomselect $fitMol $selText]

    # Set the coordinates.
    set pos [$refSel get {x y z}]
    $fitSel set {x y z} $pos
    $fitSel set beta $beta
    puts "Set [$fitSel num] atoms."

    $refSel delete
    $fitSel delete
    
    set all [atomselect $fitMol all]
    $all writepdb $outName.pdb
    $all delete

    mol delete $fitMol
    mol delete $refMol
}
exit
