#!/usr/bin/perl -w

use Getopt::Long;

use strict;
use warnings;

###################################
my $orig_pattern = "(ATOM  )(.{5}).(.{4})(.)(.{3}).(.)(.{4})(.)...(.{8})(.{8})(.{8})(.{6})(.{6})";
my $chm_pattern = "(ATOM  )(.{5}).(.{4})(.)(.{4})(.)(.{4})(.)...(.{8})(.{8})(.{8})(.{6})(.{6})";

my $orig_format = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n";
my $chm_format = "%-6s%5d %-4s%1s%4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n";

my @orig_len = (6,5,4,1,3,1,4,1,8,8,8,6,6);
my @chm_len  = (6,5,4,1,4,1,4,1,8,8,8,6,6);

###################################
# Default
my $pattern = $chm_pattern;
my $format = $chm_format;
my @len = @chm_len;
###################################

my $aoff = 0;
my $roff = 0;
my $help;
my $noaid;
my $norid;
my $debug;

my @word = (
    "Usage: renumber atom serial number and residue sequence number.\n",
    "\tpdbrenum PDB|STDIN\n",
    "Options:\n",
    "\taoff=s: offset for atom number\n",
    "\troff=s: offset for residue number\n",
    "Output:\n",
    "\tSTDOUT: new PDB\n"
);
   

sub usage {
    print STDERR @word;
    exit;
}

GetOptions(
    "aoff=s" => \$aoff,
    "roff=s" => \$roff,
    "noaid" => \$noaid,
    "norid" => \$norid,
    "help" => \$help,
    "debug" => \$debug
);

&usage if defined $help;

my $atomnum = $aoff;
my $resnum = $roff;
my $oldres = -1;
while (<>) {
    #next if /^(.{6})(.{5})(.{4})(.)(.{3})(.)(.{4})/;
    next unless /^(ATOM  )(.{5}).(.{4})(.)(.{4})(.)(.{4})(.*)/;
    #next unless /$pattern/;

    #unless ($7 == $oldres) {
    unless ($7 eq $oldres) {
	++$resnum;
	$oldres = $7;
    }
    $resnum = 0 if ($resnum == 10000);
    my $rnum = sprintf "%4d", (defined $norid) ? $7 : $resnum;
    ++$atomnum;
    my $anum = sprintf "%5d", (defined $noaid) ? $2 : $atomnum;

    print "$1$anum $3$4$5$6$rnum$8\n";
    #print "$1$anum $3$4$5$6$7$8\n";
}
#print "TER\n";

