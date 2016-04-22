set nucleotide "ADE THY CYT GUA"

foreach nt $nucleotide {

    set sel [atomselect top "resname $nt and name O4'"]
    puts [$sel num]

}


