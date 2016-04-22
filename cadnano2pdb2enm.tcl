# Based on continuity of resid, separate pdb into a list of segments
# Usage: vmd -dispdev text -e cadnano2pdb2enm.tcl -args strength cutoff input output 
# Author: Chen-Yu Li <cli56@illinois.edu> 
# 2015/5/8

## Parse arguments

# Define extrabond strength.
set strength [lindex $argv 0]

# Define cutoff distance between atoms for which extrabonds will be generated.
set cutoff [lindex $argv 1]

set input [lindex $argv 2]

set output [lindex $argv 3]


## Start doing something

set seltext "name N1 or name C2 or name N3 or name C4 or name \"C5\[M\]*\" or name C6 or name N7 or name C8 or name N9"

mol load pdb $input
set out [open $output w]


set sel [atomselect top "($seltext)"]

# Cycle through indices in sel starting from lowest, and only store atoms which are within cutoff and have indices higher than the current atom.
set indices [lsort -integer [$sel get index] ]

$sel delete

foreach i $indices {
    set iatm [atomselect top "index $i"]
    set ires [$iatm get residue]
    set icoord [lindex [ $iatm get {x y z} ] 0]

    set nearAtom [atomselect top "(not residue $ires) and ($seltext) and within $cutoff of index $i"]

    foreach cindex [lsort -integer [$nearAtom get index] ] {
        if {$cindex > $i} {
            set jatm [atomselect top "index $cindex"]
            set jcoord [lindex [$jatm get {x y z}] 0]

            set distance [veclength [vecsub $icoord $jcoord]]

            puts $out [format "bond%10d%10d%10.3g%10.3g" $i $cindex $strength $distance]
        }
        # Release memory
        catch {$jatm delete}
        catch {unset jcoord distance}

    }
    # Release memory
    $iatm delete
    catch {$nearAtom delete}
    unset ires icoord
    catch {unset cindex jcoord distance}

}

close $out

exit

