for f in slow2_*.traj; do
    awk -f convertTraj.awk $f > hall_${f}
done
