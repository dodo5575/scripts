grep WEIGHT: segmentAutocorr.log | awk '{print $2}' > tmp1.txt
rm -f tmp2.txt
for x in `cat tmp1.txt`; do
    echo ${x##*.} >> tmp2.txt
done
grep WEIGHT: segmentAutocorr.log | awk '{print $3}' > tmp1.txt
paste -d " " tmp2.txt tmp1.txt > weight_index.txt
