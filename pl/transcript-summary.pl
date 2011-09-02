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
require 'pl/sub-fasta.pl';

###############################################################################
# COMMAND LINE
###############################################################################
$| = 1; # Do not buffer output
my $VERSION = 'transcript-summary.pl 1.0';

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
            'transcript=s',
            'feature=s',
            'subcmd=s',
            'format=s',
            'in=s',
            'size=i',
            'col=i',
            'fasta=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $out;
my $outfile;
my $transcript;
my $feature;
my $subcmd = "default";
my $in;
my $size = 50;
my $col = 1;
my $fasta;
my $format = "single";

if (exists $params{format})
{
  $format = $params{format};
}

if (exists $params{in})
{
  $in = $params{in};
}

if (exists $params{size})
{
  $size = $params{size};
}

if (exists $params{col})
{
  $col = $params{col};
}

if (exists $params{fasta})
{
  $fasta = $params{fasta};
}

if (exists $params{subcmd})
{
  $subcmd = $params{subcmd};
}

if (exists $params{feature})
{
  $feature = $params{feature};
}

if (exists $params{transcript})
{
  $transcript = $params{transcript};
}

if (exists $params{out})
{
  $out = $params{out};
  open ($outfile, ">", $out) or die "cannot open > $out: $!";
}
else
{
  $outfile = *STDOUT;   
}

if ($cmd eq "getsequence")
{
  unless (exists $params{in} and 
          exists $params{col} and
          exists $params{fasta})
  {
    &printError("Command $cmd needs -in -col and -fasta options");
  }
}
elsif ($cmd eq "summary")
{
  unless (exists $params{feature} 
          and exists $params{transcript})
  {
    &printError("Command $cmd needs -feature and -transcript options");
  }
}

###############################################################################
# DATA PROCESSING
###############################################################################

if ($cmd eq "summary")
{
  my $line;

  ##############################################################
  # Read all predicted transcripts.
  my @tx;
  open TX, $transcript or die "cannot open $transcript $!";
  my $i = 0;
  while ($line = <TX>)
  {
    $i++;
    my @e = split /\s+/, $line;
    my $rec = {};
    $rec->{name}       = sprintf ("transcript%05d", $i);
    $rec->{start}      = $e[1];
    $rec->{end}        = $e[2];
    $rec->{expression} = $e[3];
    push @tx, $rec;
  }
  close TX;

  ##############################################################
  # Read all features.
  my @gx;
  open GX, $feature or die "cannot open $feature $!";
  $line = <GX>; # skip the head of the first line.
  while ($line = <GX>)
  {
    my @e = split /\s+/, $line;
    my $rec = {};
    $rec->{start}   = $e[1];
    $rec->{end}     = $e[2];
    $rec->{name}    = $e[3];
    $rec->{strand}  = $e[4];
    push @gx, $rec;
  }
  close GX;

  ##############################################################
  # Find genes in transcripts.
  for (my $i = 0; $i <= $#tx; $i++)
  {
    my $t = $tx[$i];
    my $numberFeature = 0;
    my $featureString = "F";
    for (my $j = 0; $j <= $#gx; $j++)
    {
      my $g = $gx[$j];
      # Tx:         |------------------|
      # c1:  <-->
      # c2:      <---->
      # c3:              <---->
      # c4:                         <---->
      # c5:                               <---->
      # c6:     <-------------------------->
      if ($g->{end} < $t->{start}) #c1
      {
        next;
      }
      if ($g->{start} < $t->{start} 
          and $t->{start} < $g->{end} and $g->{end} <= $t->{end}) #c2
      {
        $numberFeature++;
	my $propFeature = sprintf("%.2f", ($g->{end} - $t->{start})/($g->{end} - $g->{start}));
        $featureString .= ",$g->{name}:$g->{strand}:TS:$g->{end}:$propFeature";
      }
      if ($t->{start} <= $g->{start} and $g->{start} < $t->{end} #c3
          and $t->{start} < $g->{end} and $g->{end} < $t->{end}) 
      {
        $numberFeature++;
        $featureString .= ",$g->{name}:$g->{strand}:$g->{start}:$g->{end}:1.00";
      }
      if ($t->{start} <= $g->{start} and $g->{start} < $t->{end} #c4
          and $t->{end} < $g->{end}) 
      {
        $numberFeature++;
	my $propFeature = sprintf("%.2f", ($t->{end} - $g->{start})/($g->{end} - $g->{start}));
        $featureString .= ",$g->{name}:$g->{strand}:$g->{start}:TE:$propFeature";
      }
      if ($t->{end} < $g->{start}) #c5
      {
        next;
      }
      if ($g->{start} < $t->{start} and $t->{end} < $g->{end}) #c6
      {
        $numberFeature++;
	my $propFeature = sprintf("%.2f", ($t->{end} - $t->{start})/($g->{end} - $g->{start}));
        $featureString .= ",$g->{name}:$g->{strand}:TS:TE:$propFeature";
      }
    }
    $t->{numberFeature} = $numberFeature;
    $t->{feature} = $featureString;
  }

  ##############################################################
  # Filter transcripts.
  if ($subcmd eq "singleplus")
  {
    for (my $i = 0; $i <= $#tx; $i++)
    {
      my $t = $tx[$i];

      my @e = split /,/, $t->{feature};
      if ($#e == 1)
      {
        my @e2 = split /:/, $e[1];
	my $p = $e2[4];
	if ($p > 0.9 and $e2[1] eq '+')
	{
	}
	else
	{
          delete $tx[$i];
	}
      }
      else
      {
        delete $tx[$i];
      }
    }
  }
  elsif ($subcmd eq "unannotated")
  {
    for (my $i = 0; $i <= $#tx; $i++)
    {
      my $t = $tx[$i];

      my @e = split /,/, $t->{feature};
      if ($#e == 0)
      {
        my $l = $t->{end} - $t->{start};
	if ($t->{expression} > 3 and $l > 30)
	{
	}
	else
	{
          delete $tx[$i];
	}
      }
      else
      {
        delete $tx[$i];
      }
    }
  }

  ##############################################################
  # Print out transcripts.
  for (my $i = 0; $i <= $#tx; $i++)
  {
    my $t = $tx[$i];
    if (defined $t)
    {
      print $outfile "$t->{name}\t$t->{start}\t-$t->{end}\t$t->{expression}\t$t->{numberFeature}\t$t->{feature}\n";
    }
    # print $outfile "$t->{name}\t$t->{start}\t$t->{end}\t $t->{numberFeature}\n";
  }
}
elsif ($cmd eq "getsequence")
{
  my $seq = maFastaParse ($fasta);

  ##############################################################
  # Read all positions.
  my @pos;
  open COL, $in or die "cannot open < $in $!";
  while (<COL>)
  {
    my @e = split /\s+/;
    push @pos, $e[$col-1];
  }
  close COL;

  ##############################################################
  # Print sequences.
  foreach my $p (@pos)
  {
    my $o = $p - $size - 1; 
    my $s = substr ($seq, $o, $size + 10);
    if ($format eq "single")
    {
      print $outfile "$s\n";
    }
    else 
    {
      print $outfile ">$p\n$s\n";
    }
  }
}

if (exists $params{out})
{
  close $outfile;
}
__END__
=head1 NAME

transcript-summary - summary of transcripts

=head1 VERSION

transcript-summary 1.0

=head1 SYNOPSIS

perl pl/transcript-summary.pl summary -feature gene.pos -transcript FASTQ01.bed

perl pl/transcript-summary.pl getsequence -in 1 -col 2 -size 50 -fasta NC_004350.fna

perl pl/transcript-summary.pl getsequence -format fasta -in 1 -col 2 -size 50 -fasta NC_004350.fna

=head1 DESCRIPTION

transcript-summary will help you to summarize transcripts with associated gene
annotations.

Command:
  summary - document transcript predictions from RNA-Seq. A table is created
  with columns; transcript name, start, end, coverage, number of features, and
  feature annotations. Feature annotations include feature name, strand, feature
  start, and feature end. Feature starts could be number, TS, or TE. If
  the left-most feature of a transcript starts at or after transcript start, the
  feature start would be a genomic position. Otherwise, the feature start will
  be TS. The right-most feature of a transcript ends after the transcript end,
  then the feature end will be TE.

  subcmd singleplus - transcripts with a single gene of plus strand.
  unannotated - transcripts with no annotations.

  getsequence - gets all sequences of upstream of a given size (-size). The positions
  are given by a column (-col) of an input file (-in). The genome sequence is
  in FASTA format (-fasta).

=head1 OPTIONS

=over 8

=item B<-feature> <file>

Features of a genome. Each row represents a feature: columns include feature
name, strand, start, and end.

=item B<-transcript> <file>

Transcript file is a simple BED file.

  chr1	1	192	1	1	+	1	192	255,0,0

=item B<-out> <file>

An output file name. Standard output is used unless an output file name is
specified.

=item B<-in> <file>

An input file.

=item B<-col> <number>

Use number-th column from the input file.

=item B<-size> <number>

The size of subsequences.

=item B<-fasta> <file>

The genome sequence file in FASTA.

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
