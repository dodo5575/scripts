#name scale x y samplerate

name=$1
x=`echo "$2 * $3" | bc`
y=`echo "$2 * $4" | bc`
sam=$5

"/software/vmd-1.9-x86_64/tachyon_LINUXAMD64" -fullshade -rescale_lights 0.3 -add_skylight 1.2 -aasamples $sam ${name}.dat -format TARGA -res $x $y -o ${name}.tga


