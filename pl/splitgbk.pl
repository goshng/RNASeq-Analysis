#!/usr/bin/perl -w
use strict;

my $contigFile;
my $gbkFile = $ARGV[0];
open GBK, $gbkFile or die "cannot open < $gbkFile $!";
while (<GBK>) {
  if (/^LOCUS\s+(contig\d+)\s+/) {
    $contigFile = "temp/$1.gbk";
    open CONTIG, ">", $contigFile or die "cannot open > temp/$contigFile $!";
  }
  print CONTIG $_;
  if (/^\/\//) {
    close CONTIG;
  }
}
