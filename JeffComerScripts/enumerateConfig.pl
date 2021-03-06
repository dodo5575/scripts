#!/usr/bin/env perl
die "Usage: $0 templateFile varName start end\n" unless @ARGV == 4;

my $templateFile = @ARGV[0];
my $varName = @ARGV[1];
my $start = @ARGV[2];
my $end = @ARGV[3];

# Read the template.
open(DAT, $templateFile) || die("Could not open file $index!\n");
@data=<DAT>;
#my @indexList = split(' ', $data[0]);
close(DAT);

my $n = $end - $start + 1;

# Write the config files.
for ($i = 0; $i < $n; $i++) {
    my $val = $start + $i;
    my $outFile = $templateFile . "$val.namd";
    open(OUT, ">$outFile") || die("Could not open file $outFile!\n");

    my $defLine = "set $varName $val\n";
    print OUT $defLine;
    foreach (@data) {
	print OUT $_;
    }
}
