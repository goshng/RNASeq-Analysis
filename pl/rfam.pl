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
            'rename=s',
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
elsif ($cmd eq "srnascanner")
{
}
elsif ($cmd eq "transterm")
{
}
elsif ($cmd eq "rnaz")
{
}

################################################################################
## DATA PROCESSING
################################################################################
if ($cmd eq "cmsearch")
{
  my $i = 0;
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
      if (exists $params{rename})
      {
        $i++;
        $queryName = sprintf ("%s%05d", $params{rename}, $i);
      }
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
      my $start;
      if ($a[9] eq '+')
      {
        $start = $a[7] - 1;
        print $outfile "$a[6]\t$start\t$a[8]\t$queryName\t0\t$a[9]\n";
      }
      else
      {
        $start = $a[8] - 1;
        print $outfile "$a[6]\t$start\t$a[7]\t$queryName\t0\t$a[9]\n";
      }
    }
  }
}
elsif ($cmd eq "srnascanner")
{
  my $line;
  my $i = 0;
  while ($line = <$infile>)
  {
    $i++;
    my @a = split /\s+/, $line;
    my $start;
    my $end;
    my $name;
    my $strand;
    if ($a[0] < $a[1])
    {
      $start = $a[0] - 1;
      $end = $a[1];
      $name = sprintf ("sRNAscanner%05d", $i);
      $strand = '+';
    }
    else
    {
      $start = $a[1] - 1;
      $end = $a[0];
      $name = sprintf ("sRNAscanner%05dc", $i);
      $strand = '-';
    }
    print $outfile "chr1\t$start\t$end\t$name\t0\t$strand\n";
  }
}
elsif ($cmd eq "transterm")
{
  my $line;
  while ($line = <$infile>)
  {
    # TERM 1         2847 - 2867     + F    93  -8.3 -5.29782. We use the
    if ($line =~ /\s+TERM\s(\d+)\s+(\d+)\s*-\s*(\d+)/)
    {
      my $start = $2;
      my $end = $3;
      my $name = sprintf ("Transterm%05d",$1);
      my $strand;
      if ($start < $end)
      {
        $start = $start - 1;
        $strand = '+';
      }
      else
      {
        $start = $end - 1;
        $end = $2;
        $strand = '-';
      }
      print $outfile "chr1\t$start\t$end\t$name\t0\t$strand\n";
    }
  }
}
elsif ($cmd eq "rnaz")
{
  my $i = 0;
  my $line;
  while ($line = <$infile>)
	{
    # rnaz/chr1:1043920-1044138(-).fa: SVM RNA-class probability: 0.971500
		if ($line =~ /\:(\d+)-(\d+)\(([+-])\).+probability:\s+([\d\.]+)/)
		{
		  $i++;
		  my $start = $1;
			my $end = $2;
			my $strand = $3;
			my $score = int($4*1000);
      my $name = sprintf ("RNAz%05d",$i);
			print $outfile "chr1\t$start\t$end\t$name\t$score\t$strand\n";
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

perl rfam.pl srnascanner -in file.in

=head1 DESCRIPTION

rfam.pl parses files that are created by rfam.

=head1 OPTIONS

  command: 

  cmsearch - creates a BED file from an Infernal's cmsearch output file.
  srnascanner - creates a BED file from an sRNAscanner's output file.
  transterm - creates a BED file from transterm's output file.
	rnaz - creates a BED file from a RNAz output file.

cmsearch: An Infernal's cmsearch output file contains "Query:" line. Use the
query string for annotating a BED file. The three lines include "Hit scores:"
let me know that the following non-empty lines hits. Use the E-value to decide
the hit can be included in a line for the BED file that you create.

srnascanner: sRNAscanner finds small RNAs in intergenic regions. An output file
of sRNAscanner, either RNA_position or ng_RNA_position, contains start and end
positions. We convert them to a BED file by assuming that the chromosome is
chr1. 

transterm: transterm creates an output file. The output file contains lines of
TERM 1         2847 - 2867     + F    93  -8.3 -5.29782. We use the
positions to convert an output file from transterm to a BED file. 

rnaz: We first extract DNA sequences using BEDTools's fastaFromBed. This creates
a FASTA file with a name like chr1:1043920-1044138(-). The name is used to find
the start and end positions with strandness. After BLAST and MUSCLE with all of
the known bacterial genomes, we obtain multiple sequence alignments in ClustalW
format. We use RNAz to process those alignments. Output files of RNAz contain
a line with "SVM RNA-class probability:", which is followed by a probability
value. This value can be used to rank candidate regions of small RNAs with RIT.
We grep only those lines to create a BED file, which would produce the following
line.
rnaz/chr1:1043920-1044138(-).fa: SVM RNA-class probability: 0.971500

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** GENERAL OPTIONS *****>

=item B<-out> <file>

An output file.

=item B<***** cmsearch OPTIONS *****>

=item B<-evalue> <real>

A threshold of E-value to subset hits in an cmsearch output file.

=item B<-rename>

To rename the hits.

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
