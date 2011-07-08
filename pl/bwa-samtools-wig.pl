#!/usr/bin/perl

#===============================================================================
#   Author: Sang Chul Choi, BSCB @ Cornell University, NY
#
#   File: bwa-samtools-wig.pl
#   Date: Thu Jul  7 17:02:35 EDT 2011
#   Version: 1.0
#===============================================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

$| = 1; # Do not buffer output

my $VERSION = 'bwa-samtools-wig.pl 1.0';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);        
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'wig=s',
            'bed=s',
            'genomeLength=i',
            'readLength=i'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

bwa-samtools-wig - convert the bed file to a wig file by considering read lengths

=head1 VERSION

bwa-samtools-wig 1.0

=head1 SYNOPSIS

perl map2graph.pl [-bed file] [-wig file] 
  [-genomeLength number] [-readLength number]

=head1 DESCRIPTION

bwa-samtools-wig will help you to convert the bed file that is generated by bwa
alignment to a wig file. Two required options are the genome length and read
length.

=head1 OPTIONS

=over 8

=item B<-bed> <file>

BWA alignment procedure produces a BED file with first positions where short
reads start to align in the reference genome. A line would look like this:
----
gi|24378532|ref|NC_004350.1|	20	20	n	37	+
----
I use the values at the second and sixth columns to find the reads mapped onto
the reference genome.

=item B<-wig> <file>

A wiggle file would look like this.
----
track type=wiggle_0
fixedStep chrom=chr1 start=1 step=1 span=1
0
0
0
----

=item B<-genomeLength> <number>

Length of the reference genome

=item B<-readLength> <number>

Length of the short read sequences

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

require 'pl/sub-error.pl';

my $bed;
my $wig;
my $genomeLength;
my $readLength;

if (exists $params{bed}) 
{
  $bed = $params{bed};
}
else
{
  &printError("bed file is missing");
}

if (exists $params{wig}) 
{
  $wig = $params{wig};
}
else
{
  &printError("wig file is missing");
}
if (exists $params{genomeLength}) 
{
  $genomeLength = $params{genomeLength};
}
else
{
  &printError("genome length is missing");
}
if (exists $params{readLength}) 
{
  $readLength = $params{readLength};
}
else
{
  &printError("read length is missing");
}

$genomeLength++;
my @map = (0) x $genomeLength;
$genomeLength--;

my $lineNumber = 0;
open BED, $bed or die "Could not open $bed";
while (<BED>)
{
  $lineNumber++;
}
close BED;
my $totalNumber = $lineNumber;
my $onepercentNumber = int ($totalNumber / 100);

$lineNumber = 0;
my $processedTime = 0;
my $elapsedTime = 0;
my $startTime = time; 
open BED, $bed or die "Could not open $bed";
while (<BED>)
{
  $lineNumber++; 
  my @e = split /\s+/;
  my $pos = $e[1];
  my $strand = $e[5];
  if ($strand eq '+')
  {
    for (my $i = $pos; $i < $pos + $readLength; $i++)
    {
      if (1 <= $i and $i <= $genomeLength)
      {
        $map[$i]++;
        # print "$lineNumber: $i\n";
      }
    }
  }
  elsif ($strand eq '-')
  {
    for (my $i = $pos; $i > $pos - $readLength; $i--)
    {
      if (1 <= $i and $i <= $genomeLength)
      {
        $map[$i]++;
        # print "$lineNumber: $i\n";
      }
    }
  }
  else
  {
    die "$strand must be + or -";
  }

  if ($lineNumber % $onepercentNumber == 0)
  {
    my $endTime = time; 
    my $elapsedTime = $endTime - $startTime;
    $processedTime += $elapsedTime;
    my $remainedNumber = $totalNumber - $lineNumber;
    my $remainedTime = int(($processedTime/$lineNumber) * $remainedNumber / 60);
    print STDERR "$remainedTime min to go\r";
    $startTime = $endTime; 
  }
}
close BED;

open WIG, ">", $wig or die "Could not open $wig";
print WIG "track type=wiggle_0\n";
print WIG "fixedStep chrom=chr1 start=1 step=1 span=1\n";
for (my $i = 1; $i <= $genomeLength; $i++)
{
  print WIG "$map[$i]\n";
}
close WIG;



