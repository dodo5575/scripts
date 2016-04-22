set origin {0.0 0.0 0.0}
set basis {{1 0 0} {0 1 0} {0 0 1}}
set basis1 {{-0.2 0 0} {0 -0.2 0} {0 0 0.2}}

source renderMesh.tcl
set m [mol new]


foreach zero {0} {
    set hydrant [readQuads hydrant0.raw]

    graphics $m material Opaque
    #drawMesh hydrant $m red
    drawMeshTransform hydrant $m red $origin $basis1
}

