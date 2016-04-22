use strict;
use warnings;
use Getopt::Long;

use lib ("$ENV{HOME}/perl");
use Stats;


my @data;
while (<>) {
    unshift(@data, $_);
}

my $sum = 0;
foreach my $i (0..$#data) {
    $sum += $data[$i];
    print $sum/($i+1) . "\n";
}
