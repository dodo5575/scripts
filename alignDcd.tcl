#!/usr/local/bin/vmd
# get data to compute ellipse
# DOES ALIGNMENT BEFORE MEASURING RMSD
# to use:  vmd -dispdev text  -e alignDcd.tcl
# modified by Chen-Yu Li
# cli56@illinois.edu


# define trajectory and structure files
set psffile /home/cli56/square2plate/cadnano/square2plate-1MKCl.psf
set pdbfile /home/cli56/square2plate/cadnano/square2plate-1MKCl.pdb

set dcdfile /home/cli56/square2plate/cadnano/square2plate-1MKCl-min/output/square2plate-1MKCl-min2.dcd
set outfile /home/cli56/square2plate/cadnano/square2plate-1MKCl-min/output/movie.dcd


set selection "nucleic"


# set reference selection
mol load psf $psffile pdb $pdbfile
set selref [atomselect top $selection]
puts "number of atoms in the selection: [$selref num]"

# load trajectory and wait until loading completes 
mol load psf $psffile dcd $dcdfile
set n [ molinfo top get numframes ]
puts "  Loaded $n frames."

# set dynamic selections
set seldyn [atomselect top $selection]

# define first and last frames and stride
set firstframe 0
set lastframe $n
set stride 1
#set selref [atomselect top $selection frame firstframe]

# proceed frame by frame
for {set i $firstframe} {$i <= $lastframe} {set i [ expr $i + $stride ] } {

    
    set seldyn [atomselect top $selection frame $i]
    set selAll [atomselect top all frame $i]
    # align dynamic selection to static one
    set matrix [ measure fit $seldyn $selref ]

    $selAll move $matrix

    # measure RMSD for all heavy atoms and for protein
    
    $selAll delete
    $seldyn delete
    
    puts "frame $i" 
}

pbc wrap -all -sel "not nucleic" -compound residue -center origin

animate write dcd $outfile waitfor all

quit
