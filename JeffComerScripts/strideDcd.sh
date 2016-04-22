stride=20

for f in layer1_dopc_aq{-1,2c,1}.dcd; do
    out=${f%.*}_stride${stride}.dcd
    
    ~/bin/catdcd -o $out -stride $stride $f
done
