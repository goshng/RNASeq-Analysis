#!/usr/bin/perl

# To Do:
# 1. Check the feature of the multiple input files.

#===============================================================================
#   Author: Sang Chul Choi, BSCB @ Cornell University, NY
#
#   File: map2graph.pl
#   Date: 2011-02-02
#   Version: 0.1.0
#
#   Usage:
#      map2graph [options]
#
#      Try 'map2graph -h' for more information.
#
#   Purpose: map2graph help you to convert segemehl's map output file to a
#            graph file that can be viewed by a genome viewer such as IGB, IGV,
#            etc.
#
#   Note that I started to code this based on PRINSEQ by Robert SCHMIEDER at
#   Computational Science Research Center @ SDSU, CA as a template. Some of
#   words are his not mine, and credit should be given to him. Others are my
#   script that is possible by help from Steve Hoffmann who is the author of the
#   mapping program, segemehl.
#===============================================================================

use strict;
use warnings;

#use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Temp qw(tempfile);

$| = 1; # Do not buffer output

my $VERSION = 'map2graph 0.1.0';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);        
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'map=s',
            'out_format=i',
            'genome_length=i',
            'graph=s'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

map2graph - Conversion of map to graph

=head1 VERSION

map2graph 0.1.0

=head1 SYNOPSIS

perl map2graph.pl [-h] [-help] [-version] [-man] [-verbose] [-map input_map_file] 

=head1 DESCRIPTION

map2graph will help you to convert segemehl program's output mapping files to a
graph file that can be viewed using a genome viewer such as IGV, IGB, etc.

=head1 OPTIONS

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-man>

Print the full documentation; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-map> <file>,<file>,...

Input files in segemehl output file format that contains short read sequences
with their maps on a reference genome. The first file's name is used for the
output graph file.

=item B<-fasta> <file>

Input file in FASTA format that contains the sequence data.

=item B<-genome_length> <integer>

Length of the reference genome

=item B<***** OUTPUT OPTIONS *****>

=item B<-out_format> <integer>

To change the output format, use one of the following options. If not defined,
the output format will be of IGB graph format.

1 (IGB), 2 (IGV) or 3 (UCSC Genome Browser)

=item B<***** FILTER OPTIONS *****>

=item B<-genome-len> <integer>

Length of a reference genome. This is optional.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make map2graph better.

=head1 COPYRIGHT

Copyright (C) 2011  Sang Chul Choi

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

my $genome_length = 3000000;
my $chromosomeName = "NC_000915.1";
my (@files1);

if (exists $params{genome_length}) {
    $genome_length = $params{genome_length};
}

#Check if input file exists and check if file format is correct
if (exists $params{map}) {
    @files1 = split /,/, $params{map};
    foreach my $f (@files1)
    {
        if (-e $f) {
            #check for file format
            my $format = &checkFileFormat($f);
            unless($format eq 'map') {
                &printError('input file for -map is in '.uc($format).' format not in segemehl format');
            }
        } else {
            &printError("could not find input file \"".$f."\"");
        }
    }

} else {
    &printError("you did not specify an input map file containing the mapping of reads");
}

#check if output format is possible
my $out_igv = 0;
my $out_igb = 1;
my $out_ucsc = 0;
if (exists $params{out_format}) {
    if ($params{out_format} =~ /\D/) {
        &printError('output format option has to be an integer value');
    } 
    if ($params{out_format} >= 1 and $params{out_format} <= 3) {
        &printError('output format option has to be inclusively between 1 and 3.');
    } 
    if ($params{out_format} == 1)
    {
      $out_igv = 1;
      $out_igb = 0;
      $out_ucsc = 0;
    }
    elsif ($params{out_format} == 2)
    {
      $out_igv = 0;
      $out_igb = 1;
      $out_ucsc = 0;
    }
    elsif ($params{out_format} == 3)
    {
      $out_igv = 0;
      $out_igb = 0;
      $out_ucsc = 1;
    }
}

#
################################################################################
## DATA PROCESSING
################################################################################
#

my $filename = $files1[0];
while($filename =~ /[\w\d]+\.[\w\d]+/) {
    $filename =~ s/\.[\w\d]+$//;
}

#create filehandles for the output data
my ($fhgraph,$fhigv);
my ($filenamegraph);

$filenamegraph = $filename.'.graph';
open($fhgraph,">".$filenamegraph) or &printError('cannot open output file');

#progress bar
my $numlines = 0;
my ($progress,$counter,$part);
$progress = $counter = $part = 1;
if(exists $params{verbose}) {
    print STDERR "Estimate size of input data for status report (this might take a while for large files)\n";
    $numlines = 0;
    foreach my $f (@files1)
    {
        $numlines += &getLineNumber($f);
    }
    print STDERR "\tdone\n";

    #for progress bar
    $progress = 0;
    $counter = 1;
    $part = int($numlines/100);
}


#parse input data
print STDERR "Parse and process input data\n" if(exists $params{verbose});
my $numseqs = 0;
my $count = 0;
my (%hits);
my (@elements);
#################################################################################
# Code of mine.
#################################################################################

# Count multiple hits of a short read.
my $files1_separated_by_a_space = join " ", @files1;
open(FILE,"cat $files1_separated_by_a_space | grep \"^>\" | cut -f1 | sort | uniq -c |") or die "ERROR: Could not open one of files $files1_separated_by_a_space: $! \n";
while (<FILE>) {
   chomp();
# 0 1                2
# -----
#   1 >SRR031223.10002
#   2 >SRR031223.10010
   @elements = split /\s+/;
   $hits{$elements[2]} = $elements[1];
}
close(FILE);
#################################################
# Just a check of the %hits hash.
# -------------------------------
# while ( my ($key, $value) = each(%hits) ) {
#     print "$key => $value\n";
# }
#################################################

# Two graphs are initialized.
my @graphplus = (0) x $genome_length;
my @graphminus = (0) x $genome_length;

open(FILE,"cat $files1_separated_by_a_space | perl -p -e 's/\r/\n/g' |") or die "ERROR: Could not open one of files $files1_separated_by_a_space : $! \n";

while(<FILE>) {
   if (/^\#/) {
       next;
   }
   chomp();
# The following is the column from sege
# 0 "descr"	
# 1 "semi global alignment edist"	
# 2 "seed score"	
# 3 "seed Evalue"	
# 4 "seed qstart"	
# 5 "seed qend"	
# 6 "semi global alignment matches"	
# 7 "semi global alignment mismatches"	
# 8 "semi global alignment insertions"	
# 9 "semi global alignment deletions"	
# 10 "strand"	
# 11 "start of semi global alignment in subject sequence"	
# 12 "end of semi global alignment in subject sequence"	
# 13 "sequenceidx"
#>SRR031223.15	12	104	0.0000	20	123	116	9	2	1	+	96928	97053	>NC_000915.1	S1;M1;S2;M2;S2;M2;S1;M1;S1;M1;S1;M2;I1;M1;S1;M104;D2;M2;
#
# A possible algorithm:
# ---------------------
# Find numbers of hits for each read. Prepare two arrays of size being equal to
# the number of nucleotides of a reference genome. Reset all the elements of the
# arrays to zero. One of arrays is for the leading strand and another for the
# lagging strand. For each map of a read, locate the range of the map to increase
# the position of the array depending on strand by 1/hits where hits is the number
# of multiple hits of the read. Repeat this procedure for all the mapped reads.
# The two arrays are the result of the graph computation algorithm.

    @elements = split /\t/;
    # print $elements[0], "\t"; # "descr"	
    # print $elements[10], "\t"; # strand
    # print $elements[11], "\t"; # start
    # print $elements[12], "\n"; # end
    for (my $i = $elements[11]; $i <= $elements[12]; $i++)
    {
        my $numberOfHits = $hits{$elements[0]};
        my $unitScore = 1/$numberOfHits;
        if ($elements[10] eq '+')
        {
            $graphplus[$i] += $unitScore; 
        }
        elsif ($elements[10] eq '-')
        {
            $graphminus[$i] += $unitScore; 
        }
    }
}

print STDERR "\r\tdone          \n" if(exists $params{verbose});
close(FILE);

# print out scores along the genome.
if ($out_igv == 1)
{
    my $filenameigv = $filename.'.igv';
    open($fhigv,">".$filenameigv) or &printError('cannot open output file');
    print $fhigv "#type=CHIP\n";
    print $fhigv "Chromosome\tStart\tEnd\tFeature\tForward\tReverse\n";
    for (my $i = 0; $i < $genome_length; $i++)
    {
        my $nexti = $i + 1;
        print $fhigv "$chromosomeName\t$i\t$nexti\toperon\t", $graphplus[$i], "\t", $graphminus[$i], "\n";
    }
    close($fhigv);
}
elsif ($out_igb == 1)
{
    my ($fhigbPlus, $fhigbMinus);
    my $filenameigbPlus = $filename.'-plus.sgr';
    open($fhigbPlus,">".$filenameigbPlus) or &printError('cannot open output file');
    my $filenameigbMinus = $filename.'-minus.sgr';
    open($fhigbMinus,">".$filenameigbMinus) or &printError('cannot open output file');
    for (my $i = 0; $i < $genome_length; $i++)
    {
        my $nexti = $i + 1;
        print $fhigbPlus "$chromosomeName\t$nexti\t", $graphplus[$i], "\n";
        print $fhigbMinus "$chromosomeName\t$nexti\t", $graphminus[$i], "\n";
    }
    close($fhigbPlus);
    close($fhigbMinus);
}
elsif ($out_ucsc == 1)
{
}
####################################################################
# Just a test of the output
# -------------------------
# for (my $i = 0; $i < $genome_length; $i++)
# {
#     print $fhgraph $graphplus[$i], "\t", $graphminus[$i], "\n";
# }
####################################################################
close($fhgraph);

##
#################################################################################
### MISC FUNCTIONS
#################################################################################
##

sub printError {
    my $msg = shift;
    print STDERR "ERROR: ".$msg.".\n\nTry \'map2graph -h\' for more information.\nExit program.\n";
    exit(0);
}

sub getLineNumber {
    my $file = shift;
    my $lines = 0;
    open(FILE,"perl -p -e 's/\r/\n/g' < $file |") or die "ERROR: Could not open file $file: $! \n";
    $lines += tr/\n/\n/ while sysread(FILE, $_, 2 ** 16);
    close(FILE);
    return $lines;
}


sub checkFileFormat {
    my $file = shift;

    open(FILE,"perl -p -e 's/\r/\n/g' < $file |") or die "ERROR: Could not open file $file: $! \n";
    while (<FILE>) {
    }
    close(FILE);

    my $format = 'map';
    return $format;
}

