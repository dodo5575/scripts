#!/usr/bin/env perl
die "Usage: $0 inputFileName\n" unless @ARGV >= 1;

for ($i = 0; $i < @ARGV; $i++) {
  my $source = @ARGV[$i];
  my $dest = $source . ".force";

  print "Extracting forces from file `$source'...\n";

  die "$source does not exist!\n" unless -e $source;
  #die "$dest exists!\n" if -e $dest;

  my $awkCmd = "'{print \$3, \$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11}'";
  system "grep \"TCLFORCES \" $source | awk $awkCmd > $dest";
  print "Results written to $dest\n"
}
