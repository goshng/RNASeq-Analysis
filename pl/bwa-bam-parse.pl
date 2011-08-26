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
my $VERSION = 'bwa-bam-parse.pl 1.0';

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
            'in=s',
            'out=s',
            'reads=s',
            'at=i',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $at;
my $out;
my $outfile;
my $in;
my $infile;
my $reads;

if (exists $params{at}) 
{
  $at = $params{at};
}

if (exists $params{reads}) 
{
  $reads = $params{reads};
}

if (exists $params{out})
{
  $out = "$params{out}";
  open ($outfile, ">", $out) or die "cannot open > $out: $!";
}
else
{
  $outfile = *STDOUT;   
}

if (exists $params{in})
{
  $in = "$params{in}";
  open ($infile, "<", $in) or die "cannot open > $in: $!";
}
else
{
  $infile = *STDIN;   
}


if ($cmd eq "find")
{
  unless (exists $params{reads}) 
  {
    &printError("command find needs -reads");
  }
}
elsif ($cmd eq "list")
{
  unless (exists $params{at}) 
  {
    &printError("command list needs -at");
  }
}

###############################################################################
# DATA PROCESSING
###############################################################################
if ($cmd eq "list")
{
  # HWI-ST397_0000:3:1103:16035:134981#CGATGT     0
  # gi|24378532|ref|NC_004350.1|  1       15      15M     *       0       0
  # TTTTGTTTTTGTTTT       FBFFFDGGEGGDGGB XT:A:U  NM:i:0  X0:i:1  X1:i:7  XM:i:0
  # XO:i:0        XG:i:0  MD:Z:15
  while (<$infile>)
  {
    chomp;
    my @e = split /\t/;
    my $name = $e[0];
    my $flag = $e[1];
    my $chr = $e[2];
    my $pos = $e[3];
    my $mapQuality = $e[4];
    if ($pos == $at)
    {
      print $outfile "$_\n";
    }
    elsif ($pos < $at)
    {
      next;
    }
    elsif ($pos > $at)
    {
      last;
    }
  }
}
elsif ($cmd eq "find")
{
  my @shortReads;
  open READS, $reads or die "cannot open < $reads $!";
  while (<READS>)
  {
    chomp;
    my @e = split /\t/;
    my $name = $e[0];
    push @shortReads, $name;
  }
  close READS;

  # HWI-ST397_0000:3:1103:16035:134981#CGATGT     0
  # gi|24378532|ref|NC_004350.1|  1       15      15M     *       0       0
  # TTTTGTTTTTGTTTT       FBFFFDGGEGGDGGB XT:A:U  NM:i:0  X0:i:1  X1:i:7  XM:i:0
  # XO:i:0        XG:i:0  MD:Z:15
  my $s = 0;
  while (<$infile>)
  {
    $s++;
    chomp;
    my @e = split /\t/;
    my $name = $e[0];
    my $flag = $e[1];
    my $chr = $e[2];
    my $pos = $e[3];
    my $mapQuality = $e[4];
    if (@shortReads)
    {
      for (my $i = 0; $i <= $#shortReads; $i++)
      {
        if ($name eq $shortReads[$i])
        {
          print $outfile "$_\n";
          splice(@shortReads, $i, 1);
          last;
        }
      }
    }
    else
    {
      last;
    }
    if ($s % 100000 == 0)
    {
      print STDERR "Read $s\r";
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

bwa-bam-parse - Subsample a FASTQ file

=head1 VERSION

bwa-bam-parse 1.0

=head1 SYNOPSIS

samtools view bam.file | perl bwa-bam-parse.pl list -at position
samtools view bam.file | perl bwa-bam-parse.pl find -reads sam.file

=head1 DESCRIPTION

bwa-bam-parse will help you to parse BAM file.

=head1 OPTIONS

command: sample, cut
  list - list short reads aligned at a position
  find - find bam.file for short reads in sam.file 

=over 8

=item B<-reads> <file>

sam.file

=item B<-at> <number>

A position at which reads start

=item B<-out> <file>

An output file name.

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
