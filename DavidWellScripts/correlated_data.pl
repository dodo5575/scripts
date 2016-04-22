use strict;
use warnings;
use Getopt::Long;

use lib ("$ENV{HOME}/perl");
use Stats;

### PARAMETERS ###
#
my $n_default = 10000;
my $tau_default = 50;
my $a_default = 1e-1;


my $n;
my $tau;
GetOptions("num=i" => \$n, "tau=i" => \$tau, "a=f" => \$a);

$n = $n_default unless defined $n;
$tau = $tau_default unless defined $tau;
$a = $a_default unless defined $a;

my @data;
for (my $i = 0; $i < $n; $i++) {
    my $mean = 0;
    my $mean_curr = mean(@data);
#     for (my $j = 0; $j < $i; $j++) {
# 	$mean += $a * exp(($j-$i)/$tau) * ($data[$j] - $mean_curr);
#     }
    $mean += ($i > 0) ? $a * $data[$i-1] : 0;
    
    my $datapt = gaussian_rand() + $mean;
    push(@data, $datapt);
    
    #last if abs($mean) > 100;
}

foreach my $datapt (@data) {
    print "$datapt\n";
}


sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
    return wantarray ? ($g1, $g2) : $g1;
}
