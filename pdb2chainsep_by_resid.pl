#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

###################################
my $orig_pattern = "(ATOM  )(.{5}).(.{4})(.)(.{3}).(.)(.{4})(.)...(.{8})(.{8})(.{8})(.{6})(.{6})";
my $chm_pattern = "(ATOM  )(.{5}).(.{4})(.)(.{4})(.)(.{4})(.)...(.{8})(.{8})(.{8})(.{6})(.{6})";

my $pattern = $chm_pattern;
my $orig_format = "%-6s%5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n";
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
my $ich = 0;
my $ch;
my $old_rid;
my $fh;
while (<>) {
    next unless /^(ATOM  )(.{5}).(.{4})(.)(.{4})(.)(.{4})(.*)/;
    my $rid = $7;

    if (!defined $fh) {
	open $fh, ">dna_chain_$ich.pdb";
	$ch = chr(ord('A')+$ich%26);
	print STDERR "chain=$ch.\n";
    }
    elsif (($rid-$old_rid) != 1 && ($rid-$old_rid) != 0) {
	print $fh "TER\n";
	close $fh;
	++$ich;
	open $fh, ">dna_chain_$ich.pdb";
	$ch = chr(ord('A')+$ich%26);
	print STDERR "chain=$ch.\n";
    }

    $old_rid = $rid;
    my $rnum = sprintf "%4d", $rid;
    my $anum = sprintf "%5d", $2;

    print $fh "$1$anum $3$4$5$ch$rnum$8\n";
    #printf $fh $orig_format, $1,$2,$3,$4,$5,$ch,$7,$8,$9,$10,$11,$12,$13;
}
print $fh "TER\n";
