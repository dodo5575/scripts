#!/usr/bin/env perl
die "Usage: $0 configFile jobs machine\n" unless @ARGV == 3;

my $execFile = "/home/jcomer/bin/browntown";
my $prefix = @ARGV[0];
my $jobs = @ARGV[1];
my $machine = @ARGV[2];
my $configFile = $prefix;
my $logFile = "$prefix.log";
my $condor = "$prefix.sub";

use Cwd;
$currDir = &Cwd::cwd();

print "Running $jobs jobs on $machine.\n";

my $cmd = "cd $currDir;";
for ($i = 0; $i < $jobs; $i++) {
    $cmd .= "$execFile $configFile $configFile.$i $i > $configFile.$i.out & disown;";
}
my $sshCmd = "ssh $machine \"$cmd\"";

print "$sshCmd\n";
system($sshCmd);


