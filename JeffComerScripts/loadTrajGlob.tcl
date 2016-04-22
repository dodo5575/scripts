set nameList {neg pos nng}
set moleList {2 3 4}
set prefix "dcd/nw_layer"

foreach mole $moleList name $nameList {
    puts "MOLECULE $mole"
    mol top $mole
    animate delete all
    set dcdList [glob "${prefix}*${name}*dcd"]
   
    foreach dcd $dcdList {
	puts $dcd
	mol addfile $dcd step 10 waitfor all
    }
}
