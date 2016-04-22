for ((i=0;i<=266;i+=1)); do
    sub="s/trap_tachyon266_light.dat.tga/trap_tachyon${i}_light.dat.tga/"

    echo "sed -e $sub < trap_movie_voltage.svg > /tmp/trap.$i.svg; inkscape-export-png /tmp/trap.$i.svg; convert /tmp/trap.$i.png trap.$i.ppm"
done
