#!/usr/bin/perl -w

use Getopt::Long;
use constant PI => 4 * atan2 1, 1;

use strict;
use warnings;
use List::Util 'shuffle';

my $help;
my $namd;
my $debug;
my $cut = 10.0;
my $k=1;

my @word = (
);
   
sub usage {
    print STDERR @word;
    exit;
}

GetOptions(
    "help" => \$help,
    "namd" => \$namd,
    "k=s" => \$k,
    "cut=s" => \$cut,   # cut in angstrom
    "debug" => \$debug
);

&usage if defined $help;

my @pdb;
my $old_ch = "xxx";
my $old_resi = 9999;
my $old_resn = "xxx";
my %data;
my $nn=0;
while (<>) {
    next unless (/^(ATOM  )(.{5}).(.{4})(.)(.{4})(.)(.{4}).{4}(.{8})(.{8})(.{8})(.*)/);
    ++$nn;
    
    my ($ii,$atom,$resn,$ch,$resi,$x,$y,$z) = ($2,$3,$5,$6,$7,$8,$9,$10);
    next if ($_ =~ / H/);
    next if ($atom =~ /(P|O1P|O2P|H|')/);
    print STDERR "atom=$atom index=$nn\n";

    $atom =~ s/ //g;
    $resn =~ s/ //g;

    ## resname
    next unless ($resn =~ /(ADE|DA|THY|DT|GUA|DG|CYT|DC)/);
    #if    ($resn =~ /^(ADE|DA|A)/) { $resn = "DA"; } 
    #elsif ($resn =~ /^(THY|DT|T)/) { $resn = "DT"; }
    #elsif ($resn =~ /^(GUA|DG|G)/) { $resn = "DG"; }
    #elsif ($resn =~ /^(CYT|DC|C)/) { $resn = "DC"; }

    if ($resi != $old_resi) {
	push @pdb, {%data} if %data;

	%data = ( resn => $resn );
	$old_resi = $resi;
    }

    $data{$atom}{xyz} = [$x,$y,$z];
    $data{$atom}{index} = $nn;
}

## push last data.
push @pdb, {%data} if %data;

### calculate COM of base
for my $ri (0 .. $#pdb) {
    my @com;
    my $nnn;
    for my $ai (keys %{$pdb[$ri]}) {
	next if $ai eq "resn";
	#print "$ch..$ri..$ai\n";
	$com[0] += $pdb[$ri]{$ai}{xyz}[0];
	$com[1] += $pdb[$ri]{$ai}{xyz}[1];
	$com[2] += $pdb[$ri]{$ai}{xyz}[2];
	++$nnn;
    }
    $pdb[$ri]{com} = [$com[0]/$nnn,$com[1]/$nnn,$com[2]/$nnn];
    #printf STDERR "%f,%f,%f\n", $com[0]/$nnn,$com[1]/$nnn,$com[2]/$nnn;
}

### make network
print "[ bonds ]\n" unless defined $namd;
for my $ri (0 .. $#pdb) {
    for my $rj (($ri+1) .. $#pdb) {
	my $rd = vec_norm(vec_sub(@{$pdb[$ri]{com}},@{$pdb[$rj]{com}}));
	next if ($rd > 30);

	for my $ai (keys %{$pdb[$ri]}) {
	    next if $ai eq "resn";
	    next if $ai eq "com";
	    my @xyzi = @{$pdb[$ri]{$ai}{xyz}};
	    for my $aj (keys %{$pdb[$rj]}) {
		next if $aj eq "resn";
		next if $aj eq "com";
		my $dd = vec_norm(vec_sub(@xyzi,@{$pdb[$rj]{$aj}{xyz}}));

		next if $dd > $cut;
		if (!defined $namd) {	# gromacs
		    printf("%10d%10d%10d%10.3g%10.3g\n", 
			$pdb[$ri]{$ai}{index}, $pdb[$rj]{$aj}{index},6, $dd*0.1, $k);
		}
		else {	# namd: index is zero-based
		    printf("bond%10d%10d%10.3g%10.3g\n", 
			$pdb[$ri]{$ai}{index}-1, $pdb[$rj]{$aj}{index}-1,$k,$dd);
		}
	    }
	}
    }
}

exit;

#print "[ moleculetype ]\n";
#printf "%s\t\t%d\n","DNA",1;
#print "\n\n";
#
#print "[ atoms ]\n";
#print ";   nr       type  resnr residue  atom   cgnr     charge       mass\n";
#for my $i (0 .. $#atom) {
#    ### remove spaces
#    $atom[$i] =~ s/ //g;
#    $resn[$i] =~ s/ //g;
#
#    ### atomtype
#    my $type = $rtp{$resn[$i]}{$atom[$i]};
#    unless (defined $type) {
#	print STDERR "atom ..$atom[$i]..$resn[$i].. not found\n";
#    }
#
#    my $chg = ($atom[$i] =~ /^P/) ? -1 : 0;
#
#    # search for the element mass from the name and then
#    # store the string suitable for outputing one line for
#    # this atom type
#    my $mass;
#
#    foreach my $mass_hash ( \%masses, \%singlemasses ) {
#	$mass = find_mass $mass_hash, $atom[$i];
#	last if defined $mass;
#    }
#    die "Didn't find mass for atom type $atom[$i]\n" unless defined $mass;
#
#    printf "%5d %8s %5d %8s %8s %8d %10.3f %10.3f\n", 
#	$i+1,$type,$resi[$i],$resn[$i],$atom[$i],$i+1,$chg,$mass;
#}
#print "\n\n";
#
#print "; ENM with cut = $cut nm\n";
#print "[ bonds ]\n";
#for my $i (0 .. $#atom-1) {
#    for my $j ($i+1 .. $#atom) {
#	my $d = (defined $pbc) ? 
#	    pbc_dx($x[$i],$y[$i],$z[$i],$x[$j],$y[$j],$z[$j]) :
#	    vec_norm(vec_sub($x[$i],$y[$i],$z[$i],$x[$j],$y[$j],$z[$j]));
#
#	$d /= 10;   # nm
#	next if ($d > $cut);
#	printf "%10d%10d%10d%10.3g%10d\n", $i+1, $j+1, 1, $d, 5000;
#    }
#}


exit;



sub deg2rad { PI * $_[0] / 180.0 }

sub vec_inner {
    return ($_[0]*$_[3]+$_[1]*$_[4]+$_[2]*$_[5]);
}

sub vec_cross {
    return ($_[1]*$_[5]-$_[2]*$_[4],
            $_[2]*$_[3]-$_[0]*$_[5],
            $_[0]*$_[4]-$_[1]*$_[3]);
}

sub vec_sub {
    return ($_[0]-$_[3],$_[1]-$_[4],$_[2]-$_[5]);
}

sub vec_norm {
    return sqrt($_[0]*$_[0]+$_[1]*$_[1]+$_[2]*$_[2]);
}

sub vec_unit {
    my $len = vec_norm(@_);

    return ($_[0]/$len,$_[1]/$len,$_[2]/$len);
}
