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
my $VERSION = 'bwa-sam.pl 1.0';

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
            'mapq=i',
            'in=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $fastq;
my $out;
my $in;
my $outfile;
my $infile;
my $mapq = 0;

if (exists $params{mapq})
{
  $mapq = $params{mapq};
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

if ($cmd eq "size")
{
  my $s = 0;
  while (<$infile>)
    {
      $s++;
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
      print "$flag\t$mapQuality\t$xt\t$nm\t$x0\t$x1\t$xm\t$xo\t$xg\t$md\n";
      if ($s % 100000 == 0)
        {
          print STDERR "Reads $s\r";
        }
    }
  print STDERR "Skipped reads: $skippedRead\n"; 
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

bwa-sam - summary of BWA alignment

=head1 VERSION

bwa-sam 1.0

=head1 SYNOPSIS

perl bwa-sam.pl [command] [-in file] [-out file]

=head1 DESCRIPTION

bwa-sam will help you to summarize BWA alignment in SAM format.

Command:
  mapq  - print MAPQ values
  size  - print the number of short reads mapped
  parse - print FLAG, MAPQ, XT, X0, X1, XM, XO, XG

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

=item B<-in> <file>

If no input file is not specified, standard input is considered.

=item B<-out> <file>

An output file name. Standard output is used unless an output file name is
specified.

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
