#!/usr/bin/env perl
die "Usage: $0 start end\n" unless @ARGV == 2;

my $start = @ARGV[0];
my $end = @ARGV[1];
my $n = $end - $start + 1; 
my $bin = qdel;

# Write the config files.
for ($i = 0; $i < $n; $i++) {
    my $num = $start + $i;

    system("$bin $num");
    print("$bin $num\n");
}
