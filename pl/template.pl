#!/usr/bin/perl
###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of RNAseq Analysis.
# 
# RNAseq Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNAseq Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RNAseq Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
sub printError($);
$| = 1; # Do not buffer output
my $VERSION = 'rfam.pl 1.0';

my $cmd = ""; 
sub process {
  my ($a) = @_; 
  $cmd = $a; 
}

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help);        
GetOptions( \%params,
            'help|h',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'evalue=f',
            'in=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

################################################################################
## COMMANDLINE OPTION PROCESSING
################################################################################

my $in;
my $infile;
my $out;
my $outfile;

if (exists $params{in})
{
  $in = $params{in};
  open ($infile, "<", $in) or die "cannot open < $in: $!";
}
else
{
  $infile = *STDIN;   
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

unless (exists $params{evalue})
{
  $params{evalue} = 1e-2;
}

if ($cmd eq "")
{
  &printError("You need a command");
}
elsif ($cmd eq "cmsearch")
{
}

################################################################################
## DATA PROCESSING
################################################################################
if ($cmd eq "cmsearch")
{
  my $queryName = "";
  my $line;
  while ($line = <$infile>)
  {
    #  Query:       C0465  [CLEN=78]
    if ($line =~ /^Query:\s+(\S+)\s*\[/)
    {
      $queryName = $1;
      $line = <$infile>;
      $line = <$infile>;
      $line = <$infile>;
      last;
    }
  }

  # print $outfile "chrom\tstart\tend\tname\tscore\tstrand\n";
  while ($line = <$infile>)
  {
    last if length($line) == 0;
    my @a = split /\s+/, $line;
    last if scalar(@a) < 12;

    if ($a[3] < $params{evalue})
    {
      if ($a[9] eq '+')
      {
        print $outfile "$a[6]\t$a[7]\t$a[8]\t$queryName\t0\t$a[9]\n";
      }
      else
      {
        print $outfile "$a[6]\t$a[8]\t$a[7]\t$queryName\t0\t$a[9]\n";
      }
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

sub printError($)
{
  my $msg = shift;
  print STDERR "ERROR: ".$msg.".\n\nTry \'$VERSION -h\' for more information.\nExit program.\n";
  exit(0);
}

__END__
=head1 NAME

rfam.pl - Process rfam related files.

=head1 VERSION

rfam.pl 1.0

=head1 SYNOPSIS

perl rfam.pl cmsearch -in file.in

=head1 DESCRIPTION

rfam.pl parses files that are created by rfam.

=head1 OPTIONS

  command: 

  cmsearch - creates a BED file from an Infernal's cmsearch output file.

cmsearch: An Infernal's cmsearch output file contains "Query:" line. Use the
query string for annotating a BED file. The three lines include "Hit scores:"
let me know that the following non-empty lines hits. Use the E-value to decide
the hit can be included in a line for the BED file that you create.

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-mcl> <file>

An mcl dump file. The file contains rows of strings separated by spaces. The
number of strings decreases as you go towards the end of the file.

=item B<-bed> <file>

A BED file for genes or some of the strings in the mcl dump file.

=item B<-out> <file>

An output file.

=item B<-numgenome> <numeric>

The total number of genomes or strains in the mcl dump file. If not specified,
we infer the number using the first row of the given mcl dump file.

=item B<-orthomcl>

This allows the output file of OrthoMCL. The sequence names in the output file
from OrthoMCL consist of strain or species name and gene name separated by a
vertical line.

=item B<-orthomclgroupname>

The final output file of OrthoMCL can contain names for a group or cluster at
the first column until a colon.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make rfam.pl better.

=head1 COPYRIGHT

Copyright (C) 2012 Sang Chul Choi

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
