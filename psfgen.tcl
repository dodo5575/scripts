package require psfgen
topology top_all36_na.rtf
pdbalias atom HOH O OH2
pdbalias residue HOH TIP3
foreach suff {"" 3 5} {
  pdbalias residue DA$suff ADE
  pdbalias residue DT$suff THY
  pdbalias residue DC$suff CYT
  pdbalias residue DG$suff GUA
}

mol new test.pdb

foreach S {P000 P001 P002 P003 P004 P005 SCAF} { 
  ## write temporary pdb for segments
  set sel0 [atomselect top "segname $S"] 
  $sel0 writepdb tmp$S.pdb

  segment $S {pdb tmp$S.pdb}
  coordpdb tmp$S.pdb $S
  
  ## patch RNA to DNA

	

  

}

# set sel0 [atomselect top "segname SCAFS"] 
#  $sel0 writepdb tmpP00$S.pdb
# segment SCAF {pdb test.pdb}
# coordpdb test.pdb SCAF

regenerate angles dihedrals
guesscoord

writepdb out.pdb
writepsf out.psf
# exit

