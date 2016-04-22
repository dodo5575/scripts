#/bin/bash
# Author: Jeff Comer <jcomer2@illinois.edu>

for file in $@; do
    grep "ENERGY: " $file | awk '{print $14}' > $file.dat
done
