use strict;
use warnings;
use Getopt::Long;

use lib ("$ENV{HOME}/perl");
use Stats;


my @data;
while (<>) {
    unshift(@data, $_);
}

foreach my $i (0..$#data) {
    my @subset = @data[0..$i];
    my $mean = mean @subset;
    my $err = 1.96 * stderr @subset;
    print "$i $mean $err\n";
}
