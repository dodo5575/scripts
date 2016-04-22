use strict;
use warnings;

@ARGV != 2 and die "Usages: perl $0 <in> <factor>\n";
my ($in, $factor) = @ARGV;

$factor = eval $factor;	# so user can enter things like 1/100

open(IN, "< $in") or die;
while (<IN>) {
    if (/^#/ or /^$/) {
	print;
    } else {
	my ($c, $v, $dv) = split;
	print "$c " . ($factor * $v) . " " . ($factor * $dv) . "\n";
    }
}
close IN;
