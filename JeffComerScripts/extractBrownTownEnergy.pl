#!/usr/bin/env perl
die "Usage: $0 inputFileName\n" unless @ARGV >= 1;

for ($i = 0; $i < @ARGV; $i++) {
  my $source = @ARGV[$i];
  my $dest = $source . ".energy";
  
  print "Extracting forces from file `$source'...\n";

  die "$source does not exist!\n" unless -e $source;
  #die "$dest exists!\n" if -e $dest;

  system "grep \"E: \" $source | awk '{print \$3}' > $dest";
  print "Results written to $dest\n"
}
