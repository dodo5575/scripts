#for int in 1 2; do
for int in 0; do
    for dt in 2.5 5 10 20 40; do
	timestep=$(perl -e "print 1e-6*$dt")
	echo $timestep
	sed -e "s/TIMESTEP/$timestep/" -e "s/INTEGRATOR/$int/" < template_force.brown > force${int}_timestep${dt}.brown
    done
done
