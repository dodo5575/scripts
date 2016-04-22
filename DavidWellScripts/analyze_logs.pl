#!/usr/bin/perl

use strict;
use warnings;

use lib ('/home/rccarr2/Scripts/');
use DataOps;
use DataOps qw($n);

my $infile=shift;
my $write=shift;

die "usage: $0 infile column\n" unless ( defined $infile && defined $write);

my @names = qw ( TS BOND ANGLE DIHED IMPRP ELECT VDW BOUNDARY MISC KINETIC TOTAL TEMP TOTAL2 TOTAL3 TEMPAVG PRESSURE GPRESSURE VOLUME PRESSAVG GPRESSAVG);

my $col = -1;
for (my $i=0; $i< scalar @names; $i++) {
  if ( $write eq $names[$i] ) {
    $col=$i;
    last;
  }
}

die "Error: $write is not a valid column\n" if ( $col == -1 );

open(FH, "< $infile") or die "Can't open $infile: $!\n";

my $timestep;
while(<FH>) {
  chomp;
  
  if ( /Info: TIMESTEP/) {
    ($timestep) = $_ =~ /Info: TIMESTEP\s+(\d+)/;
  }
 
  my (@output) = $_ =~ /^ENERGY:\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)\s+($n+)/;
  
  next unless (defined $output[0]);

  my $ts = $output[0] * $timestep * 1e-6;

  print $ts, "\t", $output[$col], "\n";

} 

close(FH);
