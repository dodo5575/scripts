for f in ../namd/triplet_cgg_diff*.log; do
    fileName=${f##*/}
    name=${fileName%.*}

    cmd="grep 'SMD[[:space:]][[:space:]]' $f | awk '{print \$3,\$4,\$5}' > traj/$name.dat"
    echo $cmd
done
