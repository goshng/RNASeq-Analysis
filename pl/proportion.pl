#!/usr/bin/perl
# perl pl/proportion.pl 43 | wc -l
#
# Sample lines of a text file.
#
# 1/p = X.Y
# Z = (1/p-X)/X
# With probability Z 2*X is chosen.
# With probability 1-Z, one of 2*X-1 is chosen at random.
use strict;
use warnings;
my $m = 10000;
my $p = $ARGV[0];
$p /= 100;

my $X = int(1/$p);
my $Y = 1/$p - $X;
my $Z = $Y/$X;
for (my $i = 0; $i < $m; $i++)
{
  my $skip = 2*$X;
  my $r1 = rand();
  if ($r1 >= $Z)
  {
    my $r2 = int (rand(2*$X-1)) + 1;
    unless (1 <= $r2 and $r2 <= 2*$X-1)
    {
      die "$r2 must be between 1 and 2*$X-1 inclusively"
    }
    $skip = $r2;  
  }
  for (my $j = 0; $j < $skip - 1; $j++)
  {
    $i++;
  }
  print "$i\n";
}

__END__
use strict;
use warnings;
my $m = 10000;
#my $proportion = 2/3;
my $proportion = $ARGV[0];
$proportion /= 100;

my $X = 1;

my $Yreal;
my $Y = int(1 / $proportion);
if ($Y == 1)
{
  $Yreal = 2 - 1/$proportion; # 0.95 1/0.95=1.052632; 1:1 - 0.052632, 2:0.052632
}
else
{
  $Yreal = (1/$proportion - $Y)/$Y;
  $Y = $Y*2 - 1;
}

print STDERR "p: $proportion\n";
print STDERR "Y: $Y\n";
print STDERR "Yreal: $Yreal\n";
for (my $i = 0; $i < $m; $i++)
{
  my $skip;

  if ($Y == 1)
  {
    my $r1 = int (rand($Y-$X+1)) + $X;
    unless ($X <= $r1 and $r1 <= $Y)
    {
      die "$r1 must be between $X and $Y inclusively"
    }
    $skip = $r1;
    if ($r1 == $Y)
    {
      my $r2 = rand();
      if ($r2 < $Yreal)
      {
        # No code.
      }
      else
      {
        $skip++;
      }
    }

    for (my $j = 0; $j < $skip - 1; $j++)
    {
      $i++;
    }
    print "$i\n";
  }
  elsif ($Y > 1)
  {
    my $r1 = rand();
    if ($r1 < $Yreal)
    # if ($r1 < 0.04857143)
    {
      # no code.
      $skip = $Y + 1;
    }
    else
    {
      $skip = int (rand($Y-$X+1)) + $X;
    }
    for (my $j = 0; $j < $skip - 1; $j++)
    {
      $i++;
    }
    print "$i\n";
  }
  else
  {
    die "assert here";
  }
}
