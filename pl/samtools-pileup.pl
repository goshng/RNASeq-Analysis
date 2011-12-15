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
my $VERSION = 'samtools-pileup.pl 1.0';

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
            'refgenome=s',
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
my $refgenome;
my $genomeLength;

if (exists $params{refgenome})
{
  $refgenome = $params{refgenome};
}
else
{
  &printError("Reference genome is missing");
}

=cut
if (exists $params{genomeLength}) 
{
  $genomeLength = $params{genomeLength};
}
else
{
  &printError("genome length is missing");
}
=cut

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

###############################################################################
# DATA PROCESSING
###############################################################################

if ($cmd eq "test")
{
  my $s = 0;
  while (<$infile>)
    {
      $s++;
      chomp;
      my @e = split /\t/;

      for (my $i = 0; $i <= $#e; $i++)
        {
          print $outfile "[$i]: $e[$i]\n";
        }
      last;
    }
}

if ($cmd eq "wiggle")
{
  my %fl = peachFastaLength ($refgenome);
=cut
  foreach my $i (keys %fl)
  {
    print "$i: $fl{$i}\n";
  }
=cut
  my %maps;
  foreach my $i (keys %fl)
  {
    my $l = $fl{$i};
    my @map = (0) x ($l + 1);
    $maps{$i} = [ @map ];
  }

  my $s = 0;
  while (<$infile>)
    {
      $s++;
      chomp;
      my @e = split /\t/;
      my $chromosome = $e[0];
      my $pos = $e[1];
      my $reads = $e[3];
      $maps{$chromosome}[$pos] = $reads;

      if ($s % 100000 == 0)
        {
          print STDERR "Line $s\r";
        }
    }

  print $outfile "track type=wiggle_0\n";
  foreach my $i (keys %fl)
  {
    print $outfile "fixedStep chrom=chr1 start=1 step=1 span=1\n";
    my $l = scalar @{ $maps{$i} };
    for (my $j = 1; $j < $l; $j++)
    {
      print $outfile "$maps{$i}[$j]\n";
    }
  }
}
elsif ($cmd eq "contigwiggle")
{
  # Test this command with a pileup file:
  # output/smu86/1/bwa/FASTQ040.pileup

  my %fl = peachFastaLength ($refgenome);
=cut
  foreach my $i (keys %fl)
  {
    print "$i: $fl{$i}\n";
  }
=cut
  my %maps;
  foreach my $i (keys %fl)
  {
    my $l = $fl{$i};
    my @map = (0) x ($l + 1);
    $maps{$i} = [ @map ];
  }

  my $s = 0;
  while (<$infile>)
    {
      $s++;
      chomp;
      my @e = split /\t/;
      my $chromosome = $e[0];
      my $pos = $e[1];
      my $reads = $e[3];
      $maps{$chromosome}[$pos] = $reads;

      if ($s % 100000 == 0)
        {
          print STDERR "Line $s\r";
        }
    }

  print $outfile "track type=wiggle_0\n";
  foreach my $i (keys %fl)
  {
    print $outfile "fixedStep chrom=chr1 start=1 step=1 span=1\n";
    my $l = scalar @{ $maps{$i} };
    for (my $j = 1; $j < $l; $j++)
    {
      print $outfile "$maps{$i}[$j]\n";
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

samtools-pileup - convert pileups to a wiggle file

=head1 VERSION

samtools-pileup 1.0

=head1 SYNOPSIS

perl samtools-pileup.pl wiggle -genomeLength 1000000 -in fastq.pileup -out fastq.wig

=head1 DESCRIPTION

samtools-pileup will help you to summarize a pileup file.

Command:
  wiggle - convert the input pileup file to a wiggle file

=head1 OPTIONS

=over 8

=item B<-in> <file>

A pileup file from samtools. If no input file is not specified, standard input
is considered. 

=item B<-out> <file>

An output file name. Standard output is used unless an output file name is
specified.

=item B<-genomeLength> <number>

A length of the reference genome.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message RNASeq Analysis at codaset dot
com repository so that I can make RNASeq Analysis better.

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
