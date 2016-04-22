# Author: Jeff Comer <jcomer2@illinois.edu>
# Build a silica system of a given density

# Parameters:
set targetDensity 1.8
set membraneVolume 729856.9763
set silicaUnitCellMass 240.337
# Output:
set outName sio

set nu [expr {$targetDensity*$membraneVolume/$silicaUnitCellMass}]
set n [expr {int(ceil(sqrt($nu)))}]
package require inorganicbuilder
package require psfgen
inorganicBuilder::initMaterials
set box [inorganicBuilder::newMaterialBox SiO2 {0 0 0} [list $n $n 1]]
inorganicBuilder::buildBox $box $outName

exit
