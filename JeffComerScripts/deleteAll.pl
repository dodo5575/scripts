#!/usr/bin/env perl
die "Usage: $0 start end\n" unless @ARGV == 2;

my $bin = "qdel";
#my $prefix = @ARGV[0];
my $start = @ARGV[0];
my $end = @ARGV[1];

my $n = $end - $start + 1;

# Write the config files.
for ($i = 0; $i < $n; $i++) {
    my $val = $start + $i;
    system("$bin $val");
}
