
set outBinData0 [binary format if3f3f3 1 {-0.1 0.2 0.3} {-1 -2 -3} {10 20 30}]
set outBinData1 [binary format if3f3f3 2 {-0.5 0.6 0.7} {-1 -2 -3} {10 20 30}]
set outBinData2 [binary format if3f3f3 3 {-0.8 0.9 1.0} {-1 -2 -3} {10 20 30}]

#set outBinData0 [binary format f3f3f3 {-0.1 0.2 0.3} {-1 -2 -3} {10 20 30}]
#set outBinData1 [binary format f3f3f3 {-0.2 0.6 0.7} {-2 -2 -3} {20 20 30}]
#set outBinData2 [binary format f3f3f3 {-0.3 0.9 1.0} {-3 -2 -3} {30 20 30}]

set fp [open binfile w]
fconfigure $fp -translation binary
puts -nonewline $fp $outBinData0
puts -nonewline $fp $outBinData1
puts -nonewline $fp $outBinData2
close $fp
set fp [open binfile r]
fconfigure $fp -translation binary

set inBinData [read $fp]
close $fp
puts $inBinData
set count [binary scan $inBinData if3f3f3 step pos0 pos1 force]
puts "count: $count"
puts "scan: $step $pos0 $pos1 $force"
puts $inBinData
