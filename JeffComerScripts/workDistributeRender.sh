#!/bin/bash
# Author: Jeff Comer <jcomer2@illinois.edu>

for f in $@; do
    tachyon -fullshade -add_skylight 1.5 -rescale_lights 0.6 -res 788 680 -aasamples 2 $f -o $f.tga
done
