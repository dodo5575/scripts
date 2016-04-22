#!/bin/bash
# Author: Jeff Comer <jcomer2@illinois.edu>

cutoff=2

echo "This script will delete files *.brown.#.dcd,"
echo "where # > $cutoff, in this directory and in all directories below."
echo "Shall we proceed? (y/N)"
read -s -n 1;
if [[ ! $REPLY =~ y ]]; then echo "Not deleting"; exit -1; fi

for f in `find . -iname "*.dcd"`; do
    end=${f##*.brown.}
    num=${end%.*}

    if [[ $num -gt 2 ]]; then
	echo $f
	rm -f $f
    fi
done
