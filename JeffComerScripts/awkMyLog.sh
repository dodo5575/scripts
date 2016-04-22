#!/bin/bsh

for file in $@; do
	awk '/TCL: Timestep/ {printf("%.12g %.12g\n", 1e-6*$3,$9)}' $file > $file.dist
done
