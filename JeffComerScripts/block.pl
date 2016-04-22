use strict;
use warnings;
use Getopt::Long;

use lib ("/home/david/perl");
use Stats;

my %opts;
GetOptions(\%opts, "interval=f", "print");

$opts{int} = 0.95 unless defined $opts{int};

my @data;
while (<>) {
    push(@data, $_);
}

my $j = 0;
my $n = @data;
while ($n >= 2) {
    my $c0 = ($n-1)/$n * var(@data);
    my $x = sqrt($c0 / ($n-1));
    my $dx = $x / sqrt(2*($n-1));
    
    #print "$n\n";
    #print "$j $x $dx\n";
    print 2**$j . " $x $dx\n";
    #print "$n $x $dx\n";
    #print "$j $c0\n";
    #print "@data\n\n";
    
    my @data_blocked = block(@data);
    @data = @data_blocked;
    $n = @data;
    $j++;
}

sub block {
    my @data_in = @_;
    my @data_out;
    
    foreach my $i (0..int(@data_in/2)-1) {
	my $datapt = 0.5 * ($data_in[2*$i] + $data_in[2*$i+1]);
	push(@data_out, $datapt);
    }
    
    return @data_out;
}
