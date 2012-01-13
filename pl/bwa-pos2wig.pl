#!/usr/bin/perl
###############################################################################
# Copyright (C) 2011 Sang Chul Choi
#
# This file is part of Mauve Analysis.
# 
# Mauve Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Mauve Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Mauve Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
require 'pl/sub-error.pl';

###############################################################################
# COMMAND LINE
###############################################################################
$| = 1; # Do not buffer output
my $VERSION = 'bwa-pos2wig.pl 1.0';

my $cmd = ""; 
sub process {
  my ($a) = @_; 
  $cmd = $a; 
}

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);        
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'genomeLength=i',
            'size=i',
            'mincoverage=i',
            'in=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $out;
my $in;
my $outfile;
my $infile;
my $genomeLength;
my $size = 1;
my $mincoverage = 3;

if (exists $params{out})
{
  $out = $params{out};
  open ($outfile, ">", $out) or die "cannot open > $out: $!";
}
else
{
  $outfile = *STDOUT;   
}

if (exists $params{in})
{
  $in = $params{in};
  open ($infile, "<", $in) or die "cannot open > $in: $!";
}
else
{
  $infile = *STDIN;   
}

if (exists $params{mincoverage}) 
{
  $mincoverage = $params{mincoverage};
}

if (exists $params{size}) 
{
  $size = $params{size};
}

if (exists $params{genomeLength}) 
{
  $genomeLength = $params{genomeLength};
}

if ($cmd eq "end")
{
  unless (exists $params{genomeLength}) 
  {
    &printError("command end needs -genomeLength");
  }
}

###############################################################################
# DATA PROCESSING
###############################################################################

if ($cmd eq "perfect")
{
  my $s = 0;
  while (<$infile>)
  {
  }
  print "$s\n";
}
elsif ($cmd eq "mapq")
{
  my $s = 0;
  while (<$infile>)
    {
      $s++;
      chomp;
      my @e = split /\t/;
      print "$e[4]\n";
      if ($s % 100000 == 0)
        {
          print STDERR "Reads $s\r";
        }
    }
}
elsif ($cmd eq "parse")
{
  print "FLAG\tMAPQ\tXT\tNM\tX0\tX1\tXM\tXO\tXG\tMD\n";
  my $s = 0;
  my $skippedRead = 0;
  while (<$infile>)
    {
      $s++;
      chomp;
      my @e = split /\t/;
      my $flag = $e[1];
      my $chr = $e[2];
      my $mapQuality = $e[4];
      unless ($#e > 10)
        {
          $skippedRead++; 
          next;
        }
      $e[11] =~ /XT:A:(.+)/;  my $xt = $1;
      $e[12] =~ /NM:i:(\d+)/; my $nm = $1;
      $e[13] =~ /X0:i:(\d+)/; my $x0 = $1;
      $e[14] =~ /X1:i:(\d+)/; my $x1 = $1;
      $e[15] =~ /XM:i:(\d+)/; my $xm = $1;
      $e[16] =~ /XO:i:(\d+)/; my $xo = $1;
      $e[17] =~ /XG:i:(\d+)/; my $xg = $1;
      $e[18] =~ /MD:Z:(.+)/;  my $md = $1;
      print $outfile "$flag\t$mapQuality\t$xt\t$nm\t$x0\t$x1\t$xm\t$xo\t$xg\t$md\n";
      if ($s % 1000000 == 0)
        {
          print STDERR "Reads $s\r";
        }
    }
  print STDERR "Skipped reads: $skippedRead\n"; 
}
elsif ($cmd eq "pos")
{
  my $mapq = 30;
  print $outfile "shortReadIdentifier\tchr\tshortReadStart\tshortReadEnd\tflag\tmapQuality\txt\tnm\tx0\tx1\txm\txo\txg\tmd\n"; 
  my $s = 0;
  my $skippedRead = 0;
  while (<$infile>)
    {
      $s++;
      chomp;
      # See this thread of dicussion for finding ending positions.
      # http://sourceforge.net/mailarchive/forum.php?thread_name=4A6F0237.30806%40broadinstitute.org&forum_name=samtools-help
      my @e = split /\t/;
      my $flag = $e[1];
      my $chr = $e[2];
      my $mapQuality = $e[4];
      unless ($mapQuality > $mapq and $#e > 10)
        {
          $skippedRead++; 
          next;
        }
      $e[11] =~ /XT:A:(.+)/;  my $xt = $1;
      $e[12] =~ /NM:i:(\d+)/; my $nm = $1;
      $e[13] =~ /X0:i:(\d+)/; my $x0 = $1;
      $e[14] =~ /X1:i:(\d+)/; my $x1 = $1;
      $e[15] =~ /XM:i:(\d+)/; my $xm = $1;
      $e[16] =~ /XO:i:(\d+)/; my $xo = $1;
      $e[17] =~ /XG:i:(\d+)/; my $xg = $1;
      $e[18] =~ /MD:Z:(.+)/;  my $md = $1;
      my $shortReadIdentifier = $e[0];
      my $shortReadStart = $e[3];
      $_ = $e[5];
      my $shortReadEnd = $e[3]-1;
      s/(\d+)[NMD]/$shortReadEnd+=$1/eg;
      print $outfile "$shortReadIdentifier\t$chr\t$shortReadStart\t$shortReadEnd\t$flag\t$mapQuality\t$xt\t$nm\t$x0\t$x1\t$xm\t$xo\t$xg\t$md\n"; 
      if ($s % 1000000 == 0)
        {
          print STDERR "Reads $s\r";
        }
    }
  print STDERR "Skipped reads: $skippedRead\n"; 
}
elsif ($cmd eq "rrna")
{
  my $gff;
  die "Command rrna needs a gff file" unless length($gff) > 0;

  #####################################################################
  # Find positions of rRNA from a gff file.
  my @rRNAPosition;
  open GFF, "grep \"gbkey=rRNA\" $gff|" or die "cannot open $gff $!";
  while (<GFF>)
  {
    my @e = split /\s+/;
    my $genomerRNA = $e[3];
    my $startrRNA = $e[3];
    my $endrRNA = $e[4];
    my $inforRNA = $e[8];
    $inforRNA =~ /;locus_tag=(\w+);/;
    my $namerRNA = $1;
    # print "$startrRNA\t$endrRNA\t$namerRNA\n";
    my $rec = {};
    $rec->{chr}   = $genomerRNA;
    $rec->{name}  = $namerRNA;
    $rec->{start} = $startrRNA;
    $rec->{end}   = $endrRNA;
    $rec->{count} = 0;
    unless ($rec->{start} < $rec->{end})
    {
      die "Genes start must be less than its end";
    }

    push @rRNAPosition, $rec;
  }
  close GFF;
  # 
  #####################################################################

  #####################################################################
  # Check each read.
  my %counts;
  $counts{total} = 0;
  $counts{unique} = 0;
  $counts{repeat} = 0;
  $counts{rrnaunique} = 0;
  $counts{rrnarepeat} = 0;
  $counts{unmapped} = 0;
  my $step = 0;
  while (<$infile>)
  {
    $counts{total}++;
    $step++;
    chomp;
    my @e = split /\t/;
    my $flag = $e[1];
    my $chr = $e[2];
    my $shortReadStart = $e[3];
    my $shortReadEnd = $e[3]-1;
    my $mapQuality = $e[4];
    $_ = $e[5];
    s/(\d+)[NMD]/$shortReadEnd+=$1/eg;
    unless ($#e > 10)
    {
      $counts{unmapped}++;
      next;
    }
    $e[11] =~ /XT:A:(.+)/;  my $xt = $1;
    $e[12] =~ /NM:i:(\d+)/; my $nm = $1;
    $e[13] =~ /X0:i:(\d+)/; my $x0 = $1;
    $e[14] =~ /X1:i:(\d+)/; my $x1 = $1;
    $e[15] =~ /XM:i:(\d+)/; my $xm = $1;
    $e[16] =~ /XO:i:(\d+)/; my $xo = $1;
    $e[17] =~ /XG:i:(\d+)/; my $xg = $1;
    $e[18] =~ /MD:Z:(.+)/;  my $md = $1;
    
    my $s = {};
    $s->{name}  = $e[0];
    $s->{chr}   = $chr;
    $s->{start} = $shortReadStart;
    $s->{end}   = $shortReadEnd;
    # gene  :     |---------------|
    # case 1:                 <------> 
    # case 2:  <------> 
    # case 3:        <------>
    # case 4:  <--------------------->
    unless ($s->{start} < $s->{end})
    {
      die "short read start must be less than its end";
    }

    my $onlyone = 0;
    for (my $i = 0; $i <= $#rRNAPosition; $i++)
    {
      my $g = $rRNAPosition[$i];
      my $v = $g->{count};
      #if ($g->{chr} eq $s->{chr})
      #{
        if ($g->{start} <= $s->{start} and $s->{start} <= $g->{end}
            and $g->{start} <= $s->{end} and $s->{end} <= $g->{end})
        { 
          $v++; # case 1 or 3
          $onlyone++;
        }
      #}
      $g->{count} = $v;
    }
    unless ($onlyone == 0 or $onlyone == 1)
    {
      die "A short read must be mapped to only one rRNA";
    }

    #print "$flag\t$mapQuality\t$xt\t$nm\t$x0\t$x1\t$xm\t$xo\t$xg\t$md\n";
    if ($step % 1000000 == 0)
    {
      print STDERR "Reads $step\r";
    }

    if ($mapQuality == 0)
    {
      unless ($xt eq "R")
      {
	die "Error: $mapQuality is zero and XT is not R";
      }
      if ($x0 == 1)
      {
	die "Error: $mapQuality is zero and X0 is 1";
      }
      else
      {
	# Okay!
      }
    }
    else
    {
      unless ($xt eq "U")
      {
	die "Error: $mapQuality is zero and XT is not U";
      }
      if ($x0 == 1)
      {
	# Okay!
      }
      else
      {
	die "Error: $mapQuality is not zero and X0 is not 1";
      }
    }

    if ($onlyone == 1)
    {
      if ($mapQuality == 0)
      {
	$counts{rrnarepeat} += $onlyone;
      }
      else
      {
	$counts{rrnaunique} += $onlyone;
      }
    }
    else
    {
      if ($mapQuality == 0)
      {
	$counts{repeat}++; 
      }
      else
      {
	$counts{unique}++;
      }
    }
  }
  #
  #####################################################################

  for (my $i = 0; $i <= $#rRNAPosition; $i++)
  {
    my $g = $rRNAPosition[$i];
    # print $outfile "$g->{name}\t$g->{count}\n";
  }
  my $sumRead = $counts{unique} + $counts{repeat} 
                + $counts{rrnaunique} + $counts{rrnarepeat} + $counts{unmapped};
  print $outfile "The number of unmapped reads is $counts{unmapped}\n";
  print $outfile "The number of uniquely mapped reads is $counts{unique}\n";
  print $outfile "The number of multiply mapped reads is $counts{repeat}\n";
  print $outfile "The number of short reads mapped uniquely to rRNA is $counts{rrnaunique}\n";
  print $outfile "The number of short reads mapped multiply to rRNA is $counts{rrnarepeat}\n";
  print $outfile "The sum of the 5 groups is $sumRead\n";
  print $outfile "The total number of reads is $counts{total}\n";
}
elsif ($cmd eq "end")
{
  $genomeLength++;
  my @map1 = (0) x $genomeLength;
  my @map2 = (0) x $genomeLength;
  my @map3 = (0) x $genomeLength;
  $genomeLength--;

  #####################################################################
  # Read a pos file.
  # shortReadIdentifier chr     shortReadStart  shortReadEnd    flag
  # mapQuality  xt      nm      x0      x1      xm      xo      xg      md
  my $s = 0;
  my $line = <$infile>;
  while (<$infile>)
  {
    $s++;
    chomp;
    my @e = split /\t/;
    my $readName       = $e[0];
    my $chr            = $e[1];
    my $shortReadStart = $e[2];
    my $shortReadEnd   = $e[3];
    my $flag           = $e[4];
    my $mapQuality     = $e[5];
    my $xt             = $e[6];
    my $nm             = $e[7];
    my $x0             = $e[8];
    my $x1             = $e[9];
    my $xm             = $e[10];
    my $xo             = $e[11]; 
    my $xg             = $e[12]; 
    my $md             = $e[13];

    for (my $i = $shortReadStart; $i <= $shortReadEnd; $i++)
    {
      $map1[$i]++;
    }
    for (my $i = $shortReadStart; $i < $shortReadStart + $size; $i++)
    {
      $map2[$i]++;
    }
    for (my $i = $shortReadEnd; $i > $shortReadEnd - $size; $i--)
    {
      $map2[$i]--;
    }

    # $map2[$shortReadStart]++;
    # $map2[$shortReadEnd]++;

    if ($s % 1000000 == 0)
    {
      print STDERR "Reads $s\r";
    }
  }
  #
  #####################################################################
  for (my $i = 1; $i <= $genomeLength; $i++)
  {
    if (abs($map2[$i]) > $mincoverage)
    {
      $map3[$i] = $map2[$i]/$map1[$i];
    }
    else
    {
      $map3[$i] = 0;
    }
  }

  print $outfile "track type=wiggle_0\n";
  print $outfile "fixedStep chrom=chr1 start=1 step=1 span=1\n";
  for (my $i = 1; $i <= $genomeLength; $i++)
  {
    print $outfile "$map3[$i]\n";
  }
}
elsif ($cmd eq "first")
{
  my @map3;
  my $line = <$infile>;
  $line = <$infile>;
  while ($line = <$infile>)
  {
    chomp $line;
    push @map3, $line;
  }

  $genomeLength = $#map3 + 1;

  for (my $i = 0; $i <= $#map3; $i++)
  {
    if ($map3[$i] == 0)
    {
      # print $outfile "0\t0\n";
    }
    elsif ($map3[$i] > 0)
    {
      my $d;
      my $pos = $i + 1;
      if ($pos < $genomeLength) 
      {
        my $v = $map3[$pos];
        while ($v >= 0)
        {
          $pos++;
          if ($pos < $genomeLength) 
          {
            $v = $map3[$pos];
          }
          else
          {
            $pos--;
            last;
          }
        }
        $d = $pos - $i;
      }
      else
      {
        $d = 0;
      }
      print $outfile "$map3[$i]\t$d\n";
    }
    else
    {
      my $d;
      my $pos = $i - 1;
      if ($pos >= 0) 
      {
        my $v = $map3[$pos];
        while ($v <= 0)
        {
          $pos--;
          if ($pos >= 0) 
          {
            $v = $map3[$pos];
          }
          else
          {
            $pos++;
            last;
          }
        }
        $d = $i - $pos;
      }
      else
      {
        $d = 0;
      }
      print $outfile "$map3[$i]\t$d\n";
    }
  }
}

if (exists $params{in})
{
  close $infile;
}
if (exists $params{out})
{
  close $outfile;
}
__END__
=head1 NAME

bwa-pos2wig - summary of BWA alignment

=head1 VERSION

bwa-pos2wig 1.0

=head1 SYNOPSIS

samtools view fastq.bam | perl pl/bwa-pos2wig.pl rrna -gff a.gff > sum.rrna

perl pl/bwa-pos2wig.pl end -in FASTQ01-pos.sum -out FASTQ01-end1-5.wig

perl pl/bwa-pos2wig.pl first -in FASTQ01-end1-5.wig -out FASTQ01-end1-5.wig.first

=head1 DESCRIPTION

bwa-pos2wig will help you to summarize BWA alignment in SAM format.

Command:
  mapq  - print MAPQ values
  size  - print the number of short reads mapped
  parse - print FLAG, MAPQ, XT, X0, X1, XM, XO, XG
  pos   - print short read name, chr, start, and end positions
  rrna  - check if reads are in the set of rRNA

  end   - counts the either one base pairs
  first - find the distances to the first opposite signed value. Using the
  output file from command ``end'' for a given position whose value is non-zero
  I compute the distance to the first opposite signed value. If a given position
  has a positive value, I go to the left to find the first position whose value
  is negative. Otherwise, I go to the right to find the first position whose
  value is postive.

  1  QNAME   String
  2  FLAG    Int 
  3  RNAME   String 
  4  POS     Int
  5  MAPQ    Int 
  6  CIGAR   String 
  7  RNEXT   String 
  8  PNEXT   Int 
  9  TLEN    Int
  10 SEQ     String 
  11 QUAL    String

=head1 OPTIONS

=over 8

=item B<-mapq> <number>

Minimum mapping quality to be called short read maps.

=item B<-in> <file>

If no input file is not specified, standard input is considered.

=item B<-out> <file>

An output file name. Standard output is used unless an output file name is
specified.

=item B<-gff> <file>

The genome annotation file. Command rrna needs a gff file.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make map2graph better.

=head1 COPYRIGHT

Copyright (C) 2011  Sang Chul Choi

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
