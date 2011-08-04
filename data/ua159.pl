use strict;

my @genesUnique;
my $line;
open IN, "ua159.smu109.info.csv";
while ($line = <IN>)
{
  chomp $line;
  my $i = 0;
  while ($i != -1)
  {
    $i = index $line, "SMU109", $i;
    if ($i != -1)
    {
      my @e = split (/[\|,]/, substr ($line, $i));
      my $gene = $e[0];
      push @genesUnique, $gene;
      $i += length($gene);
      print "$gene\n";
    }
  }
}
close IN;

print "=======================================\n";
open SIG, "ua159.out.compare";
while ($line = <SIG>)
{
  chomp $line;
  my @e = split /\s+/, $line;
  my $g = $e[1];
  print "For $g\n";
  foreach my $gUnique (@genesUnique)
  {
    if ($g eq $gUnique)
    {
      print "         FOUND $g : $line\n";
    }
  }
}
close SIG;
