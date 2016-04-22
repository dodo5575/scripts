#!/usr/bin/env perl
use strict;
use lib ('/home/rccarr2/Scripts/');
use DataOps;
use DataOps qw($n);

my $name = shift;
my $window = shift;

die "usage: $0 name [Windows]\n" unless (defined $name);

unless ( defined $window ) {
  $window = '../Windows.txt';
}

my $templateFile = $name;

open(FH, $templateFile) || die("Could not open file $templateFile!\n");
my @data=<FH>;
close(DAT);

open(FH, "< ${window}") or die "Can't open window file";

while(<FH>) {
  
  chomp;
  
  my ($win, $x, $y, $z, $kx, $ky, $kz) = $_ =~ /^(\d+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)/;  
  
  next unless (defined $win);

  my $outFile = $templateFile . $win . '.namd';
  open(OUT, "> $outFile") or die("Could not open file $outFile!\n");
  
  print OUT "set pos\t$win\n";
  print OUT "set springX\t$kx\n";
  print OUT "set springZ\t$kz\n";

  foreach (@data) {
    print OUT $_;
  }
}

