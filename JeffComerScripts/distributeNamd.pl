#!/usr/bin/perl

my @machines = ("tbgl-work2","tbgl-work8","tbgl-work9","tbgl-work10");
my @loads = (6,6,6,6,2,2);
my $name = "corr2_basepair_at_neg";
my $user = "jcomer";
my $workDir = "/projects/jcomer/basepair/correlation_force/namd";

my $j = 0;
my $count = 0;
foreach $mach(@machines) {
    my $namdCmd = "cd $workDir;";
    
    for ($i = 0; $i < @loads[j]; $i++) {
	my $cmd = "namd2cvs ${name}${count}.namd > ${name}${count}.log & disown;";
	$namdCmd .= $cmd;
        print "$cmd\n";
	$count++;
    }

    my $sshCmd = "ssh -l $user $mach \"$namdCmd\"";
    #print "$sshCmd\n";
    system $sshCmd;

   $j++;
}
