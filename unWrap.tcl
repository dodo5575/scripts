



proc unWrap {structPrefix dcdPrefix outPut} {

    package require pbctools

    # Input:
    set psf $structPrefix.psf
    set dcd $dcdPrefix.dcd

    mol new $psf
    mol addfile $dcd waitfor all

    pbc unwrap -all

    animate write dcd ${outPut}.dcd waitfor all top

    mol delete top

}


if {$argc < 3} {
    puts "vmd -dispdev text -e $argv0 -args structPrefix dcdPrefix outPutPrefix"
    exit
}

set structPrefix [lindex $argv 0]
set dcdPrefix [lindex $argv 1]
set outPut [lindex $argv 2]

unWrap $structPrefix $dcdPrefix $outPut

exit

