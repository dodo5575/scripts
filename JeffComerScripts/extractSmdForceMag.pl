#!/usr/bin/env perl
die "Usage: $0 inputFileName\n" unless @ARGV >= 1;

for ($i = 0; $i < @ARGV; $i++) {
  my $source = @ARGV[$i];
  my $dest = $source . ".fmag";

  print "Extracting forces from file `$source'...\n";
  die "$source does not exist!\n" unless -e $source;

# Read the input file.
  open(INPUT, $source) || die("Could not open file $source!\n");
  @data=<INPUT>;
  @data = grep(/^SMD  /, @data);
  close(INPUT);

  
  open(OUTPUT,">$dest") || die("Cannot Open File");
  foreach(@data) {
      @lin = split(" ", $_);
      my $t = @lin[1];
      my $fx = @lin[5];
      my $fy = @lin[6];
      my $fz = @lin[7];
      my $fmag = sqrt($fx*$fx + $fy*$fy + $fz*$fz);
      
      print OUTPUT "$t $fmag\n";
  }
  close(OUTPUT);
  
#  my $awkCmd = "'{print \$2,\$8}'";
 # system "grep \"SMD  \" $source | awk $awkCmd > $dest";
  #print "Results written to $dest\n"
}
