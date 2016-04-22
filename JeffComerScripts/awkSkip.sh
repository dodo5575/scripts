#!/bin/bash

# Usage: awkSkip.sh stride file0 file1 file2

i=0
for file in $@; do
    if [[ $i == 0 ]]; then
	((i++))
	stride=$file
    else
	echo "Processing \`$file' with stride $stride."
	awk "NR % $stride == 1 {print \$0}" $file > $file.$stride.skip
    fi
done
