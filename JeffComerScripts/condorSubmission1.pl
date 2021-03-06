#!/usr/bin/env perl
die "Usage: $0 configFile jobs map0 map1 map2 map3\n" unless @ARGV == 6;

my $execFile = "/home/jcomer/bin/browntown";
my $pmfFiles = "@ARGV[2], @ARGV[3], @ARGV[4], @ARGV[5], files/extend176_chl-chl.dat, files/extend176_pot-chl.dat, files/extend176_pot-pot.dat, files/coords_no_160.txt";
my $prefix = @ARGV[0];
my $jobs = @ARGV[1];
my $configFile = $prefix;
my $logFile = "$prefix.log";
my $condor = "$prefix.sub";

use Cwd;
$currDir = &Cwd::cwd();

if (-e $condor) {
    print "$condor already exists. Overwrite? ";
    <STDIN> =~ /^[yY]\n$/ or die "\n";
}


open(SCR,">$condor") || die("Cannot Open File");
print SCR "Universe = vanilla\n";
print SCR "Executable = $execFile\n";
print SCR "Log = ${configFile}.\$(Process).log\n";
print SCR "Input = $configFile\n";
print SCR "Output = ${configFile}.\$(Process).out\n";
print SCR "Error = ${configFile}.\$(Process).err\n";
print SCR "Arguments = $configFile $configFile.\$(Process) \$(Process)\n\n";
print SCR "should_transfer_files = YES\n";
print SCR "when_to_transfer_output = ON_EXIT_OR_EVICT\n";
print SCR "transfer_input_files = $pmfFiles\n";
print SCR "Queue $jobs\n";
close(SCR);
