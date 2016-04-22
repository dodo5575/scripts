for ion in pot chl; do
    for sys in ade thy gua cyt mcyt; do
	inFile=zoom_${sys}_${ion}.png
	outFile=crop_${sys}_${ion}.png

	eval "convert $inFile -crop 674x709+230+69 $outFile"
    done
done
