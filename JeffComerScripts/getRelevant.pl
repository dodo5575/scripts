#!/usr/bin/env perl
die "Usage: $0 indexFile dataFile\n" unless @ARGV == 2;

my $indexFile = @ARGV[0];
my $dataFile = @ARGV[1];

# Read the index.
open(DAT, $indexFile) || die("Could not open file $indexFile!\n");
@data=<DAT>;
my @indexList = split(' ', $data[0]);
close(DAT);

# Write the config files.
foreach (@indexList) {
    #print("$_\n");
    system("grep \"^$_\" $dataFile");
}

