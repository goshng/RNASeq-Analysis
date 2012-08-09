#!/usr/bin/perl
###############################################################################
# Copyright (C) 2011-2012 Sang Chul Choi
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
use Bio::SeqIO;
require 'pl/sub-error.pl';
$| = 1; # Do not buffer output
my $VERSION = 'find-core-gene.pl 1.0';

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
            'numgenome=i',
            'coregenome=s',
            'mcl=s',
            'orthomclgroupname',
            'orthomcl',
            'bed=s',
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

if ($cmd eq "")
{
  &printError("You need a command");
}
elsif ($cmd eq "core")
{
  unless (exists $params{coregenome}
          and exists $params{mcl}
          and exists $params{bed})
  {
    &printError("command $cmd needs options -coregenome, -mcl and -bed");
  }
}

################################################################################
## DATA PROCESSING
################################################################################
if ($cmd eq "core")
{
  # Enumerate all of the genes in the BED file.
  my %gene;
  open BED, $params{bed} or die "cannot open $params{bed} $!";
  while (<BED>)
  {
    my @a = split /\s+/;
    $gene{$a[3]} = 0; # non-core.
  }
  close BED;

  # Count the number of groups.
  my $clusterSize = 0;
  open MCL, $params{mcl} or die "cannot open $params{mcl} $!";
  while (my $line = <MCL>)
  {
    $clusterSize++;
  }
  close MCL;

  # Create a set or a hash of the core genomes.
  my %coregenome;
  if (exists $params{coregenome})
  {
    my @a = split /\s+/, $params{coregenome};
    %coregenome = map { $_ => 1 } @a;
    $params{numgenome} = scalar(keys %coregenome);
  }

  open MCL, $params{mcl} or die "cannot open $params{mcl} $!";
  my $i = 0;
  my @group;
  my @groupsize;
  my @groupname;
  while (my $line = <MCL>)
  {
    my $isCoregeneSet = 0;
    my $l = $line;

    # We could remove the first item of the line.
    my $name; 
    if (exists $params{orthomclgroupname})
    {
      if ($l =~ /^(\S+):/)
      {
        $name = $1;
      }
      else
      {
        $name = sprintf ("%05d", $i);
      }
      $l =~ s/^\S+\s+//;
    }
    else
    {
      $name = sprintf ("%05d", $i);
    }
    push @groupname, $name;

    # Collate the strains from the line.
    # Remove some genes
    $l =~ s/SRA\_\d+//g;
    if (exists $params{orthomcl})
    {
      $l =~ s/\|\S+//g;
    }
    else
    {
      $l =~ s/\_\S+//g;
      $l =~ s/\.\S+//g;
    }
    my @a = split /\s+/, $l;
    push @groupsize, scalar(@a);

    # Print all of the strains.
    my %seen = (); my @b = grep { ! $seen{$_} ++ } @a;
    my $numberOfStrain = scalar(@b);
    if ($i == 0)
    {
      unless (exists $params{numgenome})
      {
        $params{numgenome} = $numberOfStrain;
      }
      print $outfile "========================\n";
      print $outfile "Total Strains $numberOfStrain\n";
      print $outfile "========================\n";
      print $outfile join ("\n", @b);
      print $outfile "\n";
    }
    push @group, scalar(@b);

    # 
    # Check if the set of genomes in @b is equal to the core gene set.
    # 
    if (exists $params{coregenome})
    {
      my %s2 = map { $_ => 1 } @b;
      my %s3 = map { $_ => 1 } grep { $s2{$_} } keys %coregenome;
      $isCoregeneSet = scalar(keys %s3);
    }

    foreach my $k (keys %gene)
    {
      my $existGene = 0;
      if (exists $params{orthomcl})
      {
        if ($line =~ /\|$k\s+/)
        {
          $existGene = 1;
        }
      }
      else
      {
        if ($line =~ /\s+$k\s+/)
        {
          $existGene = 1;
        }
      }
       
      if ($existGene == 1)
      {
        die "The gene $k was already assigned a group $i" unless $gene{$k} == 0;
        $gene{$k} = $isCoregeneSet;
      }
    }
    $i++;
    print STDERR "$i / $clusterSize\r";
    
  }
  close MCL;

  $i = scalar(@groupsize);
  print $outfile "========================\n";
  print $outfile "Number of genes in groups $i\n";
  print $outfile "========================\n";
  foreach my $g (@groupsize)
  {
    print $outfile $g," ";
  }
  print $outfile "\n";

     
  $i = 0; foreach my $k (keys %gene) { $i++ if $gene{$k} <  $params{numgenome}; }
  print $outfile "========================\n";
  print $outfile "Non-Core genes $i\n";
  print $outfile "========================\n";
  # Print the core genes, or gene with 1.
  $i = 0;
  foreach my $k (keys %gene)
  {
    if ($gene{$k} < $params{numgenome})
    {
      print $outfile "$k\t$gene{$k}\n";
      $i++;
    }
  }

  $i = 0; foreach my $k (keys %gene) { $i++ if $gene{$k} ==  $params{numgenome}; }
  print $outfile "========================\n";
  print $outfile "Core genes $i\n";
  print $outfile "========================\n";
  foreach my $k (keys %gene)
  {
    if ($gene{$k} == $params{numgenome})
    {
      print $outfile "$k\n";
    }
  }
  print $outfile "========================\n";
  print $outfile "The number of core genes is $i.\n";
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

find-core-gene.pl - Find core genes using a mcl result file.

=head1 VERSION

find-core-gene.pl 1.0

=head1 SYNOPSIS

perl find-core-gene.pl core -orthomcl -coregenome "SMU1 SMU2" -mcl file.tab -bed file.bed

perl find-core-gene.pl core -orthomclgroupname -coregenome "SMU1 SMU2" -mcl file.tab -bed file.bed

=head1 DESCRIPTION

find-core-gene.pl parses mcl dump files.

=head1 OPTIONS

  command: 

  core - Find core genes in a BED file using a mcl dump file.

core: We test the genes given by a BED file for their membership of core or
non-core genes. We read in the bed file to have a list of genes. A row in a mcl
dump file contains a list of genes that are associated as a group of similar
genes. Using the genes names we can tell its strains. So, we can test a gene
whether it belongs to core or non-core groups. In one pass of reading the mcl
dump file, we label each row a group and count strains, and label the genes in
the BED file the group. We should be given a string for a list of core genomes
using -coregenome option. In the script, we remove genes with SRA prefix by
assuming that they do not belong to any species.

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

=item B<-coregenome> <string>

A string for a list of core genomes.

=item B<-orthomcl>

This allows the output file of OrthoMCL. The sequence names in the output file
from OrthoMCL consist of strain or species name and gene name separated by a
vertical line.

=item B<-orthomclgroupname>

The final output file of OrthoMCL can contain names for a group or cluster at
the first column until a colon. We remove the first column.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make find-core-gene.pl better.

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
