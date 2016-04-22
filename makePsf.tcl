## usage: makePsf villin-gms.pdb  villin-noWater

puts "Usage: makePsf {inPdb outFile}"

proc makePsf {inPdb outFile} {
    ## initialize psfgen
    package require psfgen             ;# load psfgen
    resetpsf                           ;# Destroys any previous attempts
    psfcontext reset                   ;# Tosses out any previously declared topology files

    psfalias                           ;# tell psfgen what resnames to expect
    topology /home/cli56/TopPar/top_all36_na.rtf ;# tells psfgen how atoms in a residue should be bound
    topology /home/cli56/TopPar/toppar_water_ions_jejoong.str
    
    ## load inPdb
    set ID [mol new $inPdb]

    ## add atom selections to psfgen
    
    for {set i 1} {$i <= 4} {incr i} { buildDNASegmentFromSel [atomselect $ID "segname SCF${i}"] SCF }
    
    for {set i 5} {$i <= 8} {incr i} { buildDNASegmentFromSel [atomselect $ID "segname STP${i}"] STP }
    buildDNASegmentFromSel [atomselect $ID "segname MGHH"] ION 
    buildDNASegmentFromSel [atomselect $ID "resname TIP3"] WT
    #buildDNASegmentFromSel [atomselect $ID "name POT"] POT
    #buildDNASegmentFromSel [atomselect $ID "name CLA"] CLA

    ## psfgen commands
    guesscoord ;# guess the coordinates of missing atoms
    regenerate angles dihedrals ;# fixes problems with patching

    ## write psf and pdb that were added in psfgen's 
    writepdb $outFile.pdb
    writepsf $outFile.psf
    
    ## cleanup in VMD
    mol delete $ID
}


##########################################################################
#### procs used by makePsf (kindly ignore, unless you are interested) ####
##########################################################################
proc buildDNASegmentFromSel { selection {segSuffix PRO} } {
    ## add $sel to psfgen
    set tmpDir /var/tmp ;# somewhere to drop temporary files
    set ID [$selection molid] ;# current VMD mol
    

    ## psfgen works by loading a pdb for each segment (a unique 4-character identifier for a group of atoms that are usually bonded together e.g. a single component from a PDB)
    ##   these can usually be identified from the "chain" field of a pdb

    ## add each chain in $sel to psfgen as its own segment
    set segnames [lsort -unique [$selection get segname]]
    foreach segname $segnames {
	set seg ${segname}
	set sel [atomselect $ID "[$selection text] and segname $segname"]
	

	## write out temporary pdb for psfgen to read
	set tmpPdb $tmpDir/tmp.pdb
	$sel writepdb $tmpPdb
	
	## add segment to psfgen using psfgen commands
	
	if {$segSuffix == "SCF" } {
		segment $seg {
	    		first none
	    		last none
	    		pdb $tmpPdb
		}
		
		## Now we patch the RNA molecule made by psfgen to make DNA
		##     By default, psfgen makes DNA.
	
		foreach {patchSeg patchResId patchResName} [join [lsort -unique [$sel get {segname resid  resname}]]] {
			if {0} {
            			patch DEO5 $patchSeg:$patchResId 
        		} else {
				patch DEOX $patchSeg:$patchResId
			}
            	
		}

		# But we also need to make the DNA periodic
        
		set first [lindex [$sel get resid] 0]
		set last [lindex [$sel get resid] end]
		patch LKNA $seg:$last $seg:$first
		
		coordpdb $tmpPdb
	
	}

	if {$segSuffix == "STP" } {
		segment $seg {
	    		first 5TER
	    		last 3TER
	    		pdb $tmpPdb
		}
		

		## Now we patch the RNA molecule made by psfgen to make DNA
		##     By default, psfgen makes DNA.
	
		foreach {patchSeg patchResId patchResName} [join [lsort -unique [$sel get {segname resid  resname}]]] {
			if {$patchResId == 1} {
            			patch DEO5 $patchSeg:$patchResId 
        		} else {
				patch DEOX $patchSeg:$patchResId
			}
            	
		}
		coordpdb $tmpPdb
	
	}

	if {$segSuffix != "SCF" && $segSuffix != "STP"} {

		segment $seg {
	    		first none
	    		last none
	    		pdb $tmpPdb
		}
		coordpdb $tmpPdb
	}

	if {$segSuffix == "DNA"} {
		## Now we patch the RNA molecule made by psfgen to make DNA
		##     By default, psfgen makes DNA.
	
		foreach {patchSeg patchResId patchResName} [join [lsort -unique [$sel get {segname resid  resname}]]] {

	    		

	    		patch DEOX $patchSeg:$patchResId
            	

		}

		# But we also need to make the DNA periodic
        
		set first [lindex [$sel get resid] 0]
		set last [lindex [$sel get resid] end]
		patch LKNA $seg:$last $seg:$first
	}

		## That was fun and easy

    }
}    

proc psfalias {} { ##
    ## Define common aliases for psfgen (so that it knows what residue and atom names mean)

    # Here's for nucleics
    foreach suff {"" 3 5} {
    pdbalias residue DA$suff ADE
    pdbalias residue DT$suff THY
    pdbalias residue DC$suff CYT
    pdbalias residue DG$suff GUA
    }


    foreach bp { GUA CYT ADE THY URA } {
	alias atom $bp "O5\*" O5'
	alias atom $bp "C5\*" C5'
	alias atom $bp "O4\*" O4'
	alias atom $bp "C4\*" C4'
	alias atom $bp "C3\*" C3'
	alias atom $bp "O3\*" O3'
	alias atom $bp "C2\*" C2'
	alias atom $bp "O2\*" O2'
	alias atom $bp "C1\*" C1'
	alias atom $bp  OP1   O1P
	alias atom $bp  OP2   O2P
    }
    alias atom THY C7 C5M

    alias atom ILE CD1 CD
    alias atom SER HG HG1
    alias residue HIS HSE

    # Heme aliases
    alias residue HEM HEME
    alias atom HEME "N A" NA
    alias atom HEME "N B" NB
    alias atom HEME "N C" NC
    alias atom HEME "N D" ND

    # Water aliases
    alias residue HOH TIP3
    alias atom TIP3 O OH2
    alias atom TIP3 OW OH2

    # Ion aliases
    alias residue K POT
    alias atom POT K POT
    alias residue ICL CLA
    alias atom ICL CL CLA
    alias residue NA SOD
    alias atom SOD NA SOD
    alias residue CA CAL
    alias atom CAL CA CAL 
    
    alias atom ATP C1* C1'
    alias atom ATP C2* C2'
    alias atom ATP C3* C3'
    alias atom ATP C4* C4'
    alias atom ATP C5* C5'
    alias atom ATP O2* O2'
    alias atom ATP O3* O3'
    alias atom ATP O4* O4'
    alias atom ATP O5* O5'
}
