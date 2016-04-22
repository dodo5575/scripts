#!/usr/bin/env perl
die "Usage: $0 prefix start end jobStart\n" unless @ARGV == 4;

my $bin = "./runbatch";
my $prefix = @ARGV[0];
my $start = @ARGV[1];
my $end = @ARGV[2];
my $jobStart = @ARGV[3];

my $n = $end - $start + 1;

# Write the config files.
for ($i = 0; $i < $n; $i++) {
    my $val = $start + $i;
    my $file = $prefix . "$val.namd";
    my $job = $jobStart + $i;
    
    system("$bin $file 16 -w 24:00:00 -d $job");
}
