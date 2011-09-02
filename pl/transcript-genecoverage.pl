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
require 'pl/sub-pos.pl';
require 'pl/sub-fasta.pl';

###############################################################################
# COMMAND LINE
###############################################################################
$| = 1; # Do not buffer output
my $VERSION = 'transcript-genecoverage.pl 1.0';

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
            'wiggle=s',
            'feature=s',
            'coveragesize=i',
            'skipsize=i',
            'minsize=i',
            'in=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $in;
my $infile;
my $out;
my $outfile;
my $wiggle;
my $feature;
my $skipsize = 50;
my $coveragesize = 500;
my $minsize = 600;

if (exists $params{skipsize})
{
  $skipsize = $params{skipsize};
}

if (exists $params{coveragesize})
{
  $coveragesize = $params{coveragesize};
}

if (exists $params{minsize})
{
  $minsize = $params{minsize};
}

if (exists $params{in})
{
  $in = $params{in};
}

if (exists $params{feature})
{
  $feature = $params{feature};
}

if (exists $params{wiggle})
{
  $wiggle = $params{wiggle};
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

if ($cmd eq "coverage" 
    or $cmd eq "mean")
{
  unless (exists $params{feature} and 
          exists $params{wiggle})
  {
    &printError("Command $cmd needs -wiggle and -feature options");
  }
}

###############################################################################
# DATA PROCESSING
###############################################################################

##############################################################
# Read the gene file.
my @genes = rnaseqPosParse ($feature); 

##############################################################
# Read the wiggle file.
my @wx;
open WIGGLE, $wiggle or die "cannot open $wiggle $!";
my $line = <WIGGLE>;
$line = <WIGGLE>;
while ($line = <WIGGLE>)
{
  chomp $line;
  push @wx, $line;
}
close WIGGLE;

if ($cmd eq "coverage")
{
  ##############################################################
  # Find genes and their coverages.
  for (my $i = 0; $i <= $#genes; $i++)
  {
    my $g = $genes[$i];
    next if ($g->{end} - $g->{start} + 1 < $minsize);
    my $start;
    my $end;
    if ($g->{strand} eq '+')
    {
      $start = $g->{start} + $skipsize - 1;
      $end = $start + $coveragesize - 1;
    }
    elsif ($g->{strand} eq '-')
    {
      $end = $g->{end} - $skipsize - 1;
      $start = $end - $coveragesize + 1;
    }
    else
    {
      die "Strand must be + or -";
    }
    unless (0 <= $start and $start <= $#wx and
            0 <= $end and $end <= $#wx)
    {
      die "$start and $end are not within the valid genomic positions";
    }

    print $outfile "$wx[$start]";
    for (my $pos = $start + 1; $pos <= $end; $pos++)
    {
      print $outfile "\t$wx[$pos]";
    }
    print $outfile "\n";
  }
}
elsif ($cmd eq "mean")
{
  ##############################################################
  # Find genes and their coverages.
  for (my $i = 0; $i <= $#genes; $i++)
  {
    my $g = $genes[$i];
    my $start = $g->{start} - 1;
    my $end = $g->{end} - 1;
    unless (0 <= $start and $start <= $#wx and
            0 <= $end and $end <= $#wx)
    {
      die "$start and $end are not within the valid genomic positions";
    }

    my $v = 0;
    for (my $pos = $start; $pos <= $end; $pos++)
    {
      $v += $wx[$pos];
    }
    $v /= ($end - $start + 1);
    print $outfile "$v\n";
  }
}


if (exists $params{out})
{
  close $outfile;
}
__END__
=head1 NAME

transcript-genecoverage - gene coverages

=head1 VERSION

transcript-genecoverage 1.0

=head1 SYNOPSIS

perl pl/transcript-genecoverage.pl coverage -feature gene.pos -transcript FASTQ01.wig

perl pl/transcript-genecoverage.pl mean -feature gene.pos -transcript FASTQ01.wig

=head1 DESCRIPTION

transcript-genecoverage will help you to extract gene coverage. Two files are
required: a wiggle file from a pileup file, and gene features.

Command:
  coverage - coverage values for genes are created. Default start position of
  the coverage is 50 downstream from the 5' end. The first parts of a fixed size
  of the coverage values are returned. Genes with minimum length are considered.

  mean - mean coverage values for genes are created.

=head1 OPTIONS

=over 8

=item B<-feature> <file>

Features of a genome. 

=item B<-wiggle> <file>

The wiggle file of coverage values.

=item B<-out> <file>

An output file name. Standard output is used unless an output file name is
specified.

=item B<-coveragesize> <number>

The size of coverage values.

=item B<-skipsize> <number>

The size of the first bases to skip.

=item B<-minsize> <number>

The minimum size of genes. Skip genes whose lengths are less than minsize.

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
