#!/usr/bin/env perl
die "Usage: $0 dir inputFileName\n" unless @ARGV >= 2;

my $dir = @ARGV[0];
for ($i = 1; $i < @ARGV; $i++) {
  my $source = @ARGV[$i];

  die "$source does not exist!\n" unless -e $source;
  my $cmd = "cd $dir; put $source";
  system "mssftp \"$cmd\"";
  #print "mssftp \"$cmd\"\n";
}
