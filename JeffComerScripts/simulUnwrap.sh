vmd -dispdev text -e analUnwrap.tcl -args layer_dopc_aq2c nw_anneal_dopc_q2c . 20 nw_layer{0,1,2}_dopc_aq2c.dcd 
vmd -dispdev text -e analUnwrap.tcl -args layer_dopc_q2c nw_anneal_dopc_q2c . 20 nw_layer{0,1,2}_dopc_q2c.dcd 
vmd -dispdev text -e analUnwrap.tcl -args layer_dopc_aq1 nw_anneal_dopc_q1 . 20 nw_layer{0,1,2}_dopc_aq1.dcd 
vmd -dispdev text -e analUnwrap.tcl -args layer_dopc_q1 nw_anneal_dopc_q1 . 20 nw_layer{0,1,2}_dopc_q1.dcd 
vmd -dispdev text -e analUnwrap.tcl -args layer_dopc_aq-1 nw_anneal_dopc_q-1 . 20 nw_layer{0,1,2}_dopc_aq-1.dcd 
vmd -dispdev text -e analUnwrap.tcl -args layer_dopc_q-1 nw_anneal_dopc_q-1 . 20 nw_layer{0,1,2}_dopc_q-1.dcd 
vmd -dispdev text -e analUnwrap.tcl -args layer_dopc_zer nw_extra_layer_zer . 20 nw_layer{0,1,2,3,4,5,6}_dopc_zer.dcd 
vmd -dispdev text -e analUnwrap.tcl -args layer_dopc_q0 nw_anneal_dopc_q0 . 20 nw_layer{0,1}_dopc_q0.dcd 
