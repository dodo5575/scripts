use strict;
use warnings;
use Getopt::Long;

use lib ("$ENV{HOME}/perl");
use Stats;

my $n; my $corr_length; my $N; my $int;
my $print;
GetOptions("sample_size=i" => \$n, "corr=i" => \$corr_length, "num=i" => \$N, "interval=f" => \$int, "print" => \$print);
(defined($n) and defined($corr_length)) and die;

$N = 10000 unless defined $N;
$int = 0.95 unless defined $int;

my @data;
while (<>) {
    push(@data, $_);
}

if (defined $n) {
    # do nothing
    next;
} elsif (defined $corr_length) {
    $n = int(@data / $corr_length);
} else {
    $n = @data;
}

($n < 1 or $N < 1) and die;

unless ($print) {
    print "n = $n  N = $N  int = $int\n";
}

my $mean_norm = mean @data;
my $stderr_norm = stderr @data;
my $err_norm = 1.96 * $stderr_norm;

my @data_out;
for (my $i = 0; $i < $N; $i++) {
    my @subdata;
    for (my $j = 0; $j < $n; $j++) {
	my $ind = int(rand $#data);
	push(@subdata, $data[$ind]);
    }
    push(@data_out, mean @subdata);
}

my $mean_bs = mean @data_out;

my @data_out_sort = sort {$a <=> $b} @data_out;
my $err_bs_up = @data_out_sort[int($int * @data_out_sort)] - $mean_bs;
my $err_bs_dn = $mean_bs - @data_out_sort[int((1 - $int) * @data_out_sort)];
my $err_bs_sym = 1.96 * stddev @data_out;

if ($print) {
    $" = "\n";
    print "@data_out\n";
} else {
    print "normal: $mean_norm +/- $err_norm\n";
    print "bootstrap: $mean_bs +/- $err_bs_up/$err_bs_dn\n";
    print "bootstrap (sym): $mean_bs +/- $err_bs_sym\n";
}
