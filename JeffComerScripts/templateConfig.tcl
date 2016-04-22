# Parameters:
set var name
set template template_eq0.namd
set outPrefix triplet_
set outSuffix _eq0.namd
set outSet {cgg gcc tac gta}

#set in [open $template]
#set templateData [read $in]
#close $in

foreach sys $outSet {
    set outFile ${outPrefix}${sys}${outSuffix}
    set out [open $outFile w]
    
    puts $out "set $var $sys"
    #puts $out $templateData
    puts $out "source $template"
        
    close $out
}
