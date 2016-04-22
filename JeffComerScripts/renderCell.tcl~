source renderMesh.tcl
set m [mol new]


foreach zero {0} {
    set elec [readQuads cell8_sub.raw]
    #set elec [readQuads cell7_electronics.raw]
    set flag [readQuads cell8_flagellum.raw]
    set nuc [readQuads cell7_nucleus.raw]
    set wall [readQuads cell7_wall.raw]
    set trod [readQuads cell8_electrodes.raw]

    graphics $m material Opaque
    drawMesh elec $m silver
    drawMesh trod $m green
    graphics $m material Glossy
    drawMesh flag $m blue2
    drawMesh nuc $m yellow
    drawMesh wall $m red
}




