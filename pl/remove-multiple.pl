# perl remove-multiple.pl < SMU86.gf3 > SMU86.gf3.new
my $chr = "x";
my $prevChr = "y";
while (<>)
{
  if (/(.+)\tGenBank\tregion/)
  {
    $prevChr = $chr;
    $chr = $1;
    unless ($chr eq $prevChr)
    {
      print;
    }
  }
  else
  {
    print;
  }
}
