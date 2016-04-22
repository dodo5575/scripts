set prefix excludeMembrane25
set inFile $prefix.dx

source vector.tcl
source gridForce.tcl

readDx g $inFile
puts "origin: $g(origin)"
puts "size: [llength $g(data)]"

ghostGrid g gg
puts "ghost origin: $gg(origin)"
puts "ghost data length: [llength $gg(data)]"
writeDx gg ${prefix}_ghost.dx

#sampleForce g ng 30 30 30 2
gridForce gg gx gy gz 
puts "force origin: $gx(origin)"
puts "force data length: [llength $gx(data)]"
writeDx gx ${prefix}_fx.dx
writeDx gy ${prefix}_fy.dx
writeDx gz ${prefix}_fz.dx

puts "Computing force magnitude!"
gridMagnitude gx gy gz gm
writeDx gm ${prefix}_fm.dx
exit



