gridScaleShift ../plate/plate_10-10-30.dx 0 1 one.dx
gridExp ../plate/plate_10-10-30.dx 1 -1 plate_exp.dx
mean=$(gridAverageMask plate_exp.dx one.dx)

size=$(gridGetSize one.dx | grep "size" | awk '{print $2}')
vol=$(gridGetSize one.dx | grep "vol" | awk '{print $2}')
echo "mean $mean size $size vol $vol"

conc0=$(perl -e "print 1.0/($vol*$mean);")
echo "$conc0 particles/AA^3"
echo "Convert to mol/l:"
units "$conc0 particles/AA^3" 'mol/l'
