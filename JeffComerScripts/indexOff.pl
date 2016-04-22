#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 indexFile directory outputFile\n" unless @ARGV == 3;
my ($indexFile, $inDir, $dest) = @ARGV;
my $timeFactor = 2e-6;
my $minEqTim = 10.0;

# Read the input file.
open(INPUT, $indexFile) || die("Could not open file `$indexFile'!\n");
open(OUTPUT,">$dest") || die("Could not open file `$dest'!");

while(<INPUT>) {
    # Continue if a comment.
    next if /^#/;
    my @lin = split;

    # Get the fields.
    my $name = $lin[0];
    my $pore = $lin[1];
    my $init = $lin[2];
    my $voltage = $lin[3];

    # Determine the input file name.
    my $inFile = "${inDir}/${name}.log";
    #print "$inFile\n";

    # Extract the off time.
    open(IN, $inFile) || die("Could not open file `$inFile'!\n");
    my @data = <IN>;
    my @off = grep(/TCLFORCESOFF/, @data);
    my @rest = grep(/ENERGY:/, @data);
    close(IN);

    # Get the length of simulation.
    my $totalTim = 0;
    if (@rest > 0) {
	my @stepLin = split(" ", $rest[-1]);
	$totalTim = $stepLin[1] * $timeFactor;
    }

    my $out;
    my $tim = -1;
    my $eqTim = 0;
    # Get the off time time.
    if (@off > 0) {
 	my @timLin = split(" ", $off[0]);
	$tim = $timLin[2] * $timeFactor;
	$eqTim = $totalTim - $tim;
    }    
    $out = "$name $pore $init $voltage $tim $eqTim\n";

    # Write the results.
    print $out;
    if ($eqTim > $minEqTim) {
	print OUTPUT $out;
    }
}

close(INPUT);
close(OUTPUT);
