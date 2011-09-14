# Author: Sang Chul Choi
# Date  : Thu Sep  8 10:17:06 EDT 2011
#
# DESCRIPTION:
# Preparation of alignments for an input file to RNAz.
# I aligned 14 genome sequences by using progressiveMauve. After preparing
# genome sequences in proper locations specified in full_alignment.xmfa I ran
# xmfa2maf to convert an XMFA formatted file to an MAF file. I prepared a set of
# integenic regions in BED format. A PHAST tool, maf_parse, used the BED and the
# MAF file to create alignments of the intergenic regions. I filtered the MAF
# file with the requirements: a) The 1st sequence must be NC_004350.gbk, b) The
# first character of each line must start with either `a' or `s', c) At least
# two sequences must be in the alignment.
#
# SYNOPSIS:
# The current directory was /Users/goshng/Documents/Projects/rnaseq/output/smutans12/1/run-mauve/output
#
# xmfa2maf full_alignment.xmfa full_alignment.maf
# maf_parse --features ../../data/feature-genome.out-intergeniconly full_alignment.maf > full_alignment.igr
# perl rnazParseMaf.pl < full_alignment.igr > full_alignment.rnaz
my @block;
my $line;
while ($line = <>)
{
  if ($line =~ /^a/)
  {
    push @block, $line;
  }
  elsif ($line =~ /^s\s+/)
  {
    push @block, $line;
  }
  elsif ($line =~ /^$/)
  {
    push @block, $line;
    if ($block[1] =~ /NC\_004350\.gbk/ and $#block > 2)
    {
      for (my $i = 0; $i <= $#block; $i++)
      {
        print $block[$i];
      }
    }
    @block = ();
  }
}
