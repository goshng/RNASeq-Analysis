
sub printError {
  my $msg = shift;
  print STDERR "ERROR: ".$msg.".\n\nTry \'$VERSION -h\' for more information.\nExit program.\n";
  exit(0);
}

1;
