while (<>) {
  if (/^seq(\d+)(.+)/) {
    my $n = int($1);
    print "chr$n$2\n"; 
  } else {
    print;
  }
}

