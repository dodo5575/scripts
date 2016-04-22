#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

###################################
my $orig_pattern = "(ATOM  )(.{5}).(.{4})(.)(.{3}).(.)(.{4})(.)...(.{8})(.{8})(.{8})(.{6})(.{6})";
my $chm_pattern = "(ATOM  )(.{5}).(.{4})(.)(.{4})(.)(.{4})(.)...(.{8})(.{8})(.{8})(.{6})(.{6})";

my $pattern = $chm_pattern;
###################################

my $col;
my $original;
my $sep = '';
my $format = "%10.5f";
my $help;

GetOptions(
    "help" => \$help
);

$pattern = $orig_pattern if defined $original;
&usage if (defined $help);

###################################
# Read Original PDB
my $index;
while (<>) {
    next unless /$orig_pattern/;
    my $atom = $3;
    ++$index;

    next unless $atom =~ "MG";
    ## extrabonds index zero-based.
    printf "bond %8d %8d %8.3f %8.3f\n", $index-1, $index, 5000, 1.94;
    printf "bond %8d %8d %8.3f %8.3f\n", $index-1, $index+3, 5000, 1.94;
    printf "bond %8d %8d %8.3f %8.3f\n", $index-1, $index+6, 5000, 1.94;
    printf "bond %8d %8d %8.3f %8.3f\n", $index-1, $index+9, 5000, 1.94;
    printf "bond %8d %8d %8.3f %8.3f\n", $index-1, $index+12, 5000, 1.94;
    printf "bond %8d %8d %8.3f %8.3f\n", $index-1, $index+15, 5000, 1.94;
}

