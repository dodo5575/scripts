#!/usr/bin/env perl
die "Usage: $0 templateFile\n" unless @ARGV == 1;

my $templateFile = @ARGV[0];
my @nucList = (ade,thy,gua,cyt);
my @ionList = (pot,chl);

# Read the template.
open(DAT, $templateFile) || die("Could not open file $index!\n");
@data=<DAT>;
#my @indexList = split(' ', $data[0]);
close(DAT);

# Write the config files.
foreach $nuc (@nucList) {
    foreach $ion (@ionList) {
	my $outFile = $templateFile . "_" . $nuc . "_" . $ion;
	open(OUT, ">$outFile") || die("Could not open file $outFile!\n");

	print OUT "set nuc $nuc\n";
	print OUT "set ion $ion\n";
	    
	foreach (@data) {
	    print OUT $_;
	}
    }
}
