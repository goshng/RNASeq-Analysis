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
require 'pl/sub-bed.pl';
require 'pl/sub-fasta.pl';
require 'pl/sub-ptt.pl';
$| = 1; # Do not buffer output
my $VERSION = 'geneontology.pl 1.0';

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
            'gff=s',
            'obo=s',
            'goa=s',
            'blast=s',
            'gene2go=s',
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
my $obo;
my $goa;
my $gene2go;
my $blast;
my $gff;

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

if ($cmd eq "gene2go")
{
  unless (exists $params{gff} and exists $params{blast}
          and exists $params{goa})
  {
    &printError("command gene2go needs options -gff, -blast, and -goa");
  }
}
elsif ($cmd eq "go2ngene")
{
  unless (exists $params{gene2go} and exists $params{obo})
  {
    &printError("command gene2go needs options -gene2go, and -obo");
  }
}

################################################################################
## DATA PROCESSING
################################################################################
if ($cmd eq "gene2go")
{
  # Parse gff file.
  my %locus2protein;
  my $geneid;
  my $locustagid;
  open GFF, $params{gff} or die "cannot open $params{gff} $!";
  while (<GFF>)
  {
    chomp;
    next if /^#/;
# NC_004350.1	RefSeq	CDS	194	1549	.	+	0	ID=NC_004350.1:dnaA:unknown_transcript_1;Parent=NC_004350.1:dnaA;locus_tag=SMU.01;note=binds to the dnaA-box as an ATP-bound complex at the origin of replication during the initiation of chromosomal replication%3B can also affect transcription of multiple genes including itself.;transl_table=11;product=chromosomal replication initiation protein;protein_id=NP_720488.1;db_xref=GI:24378533;db_xref=GeneID:1029426;exon_number=1
    # if (/locus\_tag=(\w+);.+;protein\_id=(\w+);/)

    # The S. mutans version 2 does not work with the following line.
    # if (/RefSeq\s+CDS\s+.+locus\_tag=([\w.]+);.+;protein\_id=(\w+)/)

    if (/RefSeq\s+gene\s+.+ID=(\w+);.+;locus\_tag=([\w.]+)/)
    {
      $geneid = $1; 
      $locustagid = $2;
    }
    elsif (/RefSeq\s+CDS\s+.+Parent=(\w+);.+;protein\_id=(\w+)/)
    {
      if ($geneid eq $1)
      {
        if (exists $locus2protein{$locustagid})
        {
          die "GFF: $locustagid was already added.";
        }
        $locus2protein{$locustagid} = $2;
      }
    }
  }
  close GFF;

  # Test
=comment
  foreach my $locusid (keys %locus2protein)
  {
    my $protein = $locus2protein{$locusid};
    print "$locusid\t$protein\n";
  }
  exit; 
=cut

  # Parse BLAST file.
  my %protein2uniref;
  open BLAST, $params{blast} or die "cannot open < $params{blast} $!";
  while (<BLAST>)
  {
    chomp;
    my $proteinid;
    my $uniref;
    # gi|24378533|ref|NP_720488.1|    UniRef90_Q8DWN9 100.00  452     0       0 1       452     1       452     0.0      922
    my @e = split /\t/;
    if ($e[0] =~ /(NP\w+)/)
    {
      $proteinid = $1;
    }
    elsif ($e[0] =~ /(YP\w+)/)
    {
      $proteinid = $1;
    }
    else
    {
      die "Protein ID must prefix with NP.";
    }
    if ($e[1] =~ /UniRef90\_(\w+)/)
    {
      $uniref = $1;
    }
    else
    {
      die "UniRef must prefix UniRef90.";
    }
    my $evalue = $e[10];
    if (exists $protein2uniref{$proteinid})
    {
      next;
    }
    else
    {
      my $rec = {};
      $protein2uniref{$proteinid} = $rec;
    }
    my $rec = $protein2uniref{$proteinid};
    if (exists $rec->{$uniref})
    {
      next;
      die "BLAST: $proteinid - $uniref was already added.";
    }
    $rec->{$uniref} = $evalue;
  }
  close BLAST;

=comment
  print "========== BLAST =========\n";
  foreach my $protein (keys %protein2uniref)
  {
    print $protein;
    foreach my $uniref (keys %{ $protein2uniref{$protein} })
    {
      print " $uniref:$protein2uniref{$protein}{$uniref}\t";
    }
    print "\n";
  }
  exit;
=cut

  # Parse the goa file.
  my $count = 0;
  my %uniref2goterm;
  open GOA, $params{goa} or die "cannot open < $params{goa} $!";
  while (<GOA>)
  {
    next if /^!/;
    $count++;
    chomp;
    # UniProtKB       A0A098  GSM1            GO:0005515      PMID:18510927   IPI UniProtKB:Q9XEK5        F     
    my @e = split /\s+/;
    my $goterm = $e[3];
    my $uniref = $e[1];
    if (exists $uniref2goterm{$uniref})
    {
      my $gotermArrayRef = $uniref2goterm{$uniref};
      push @{ $gotermArrayRef }, $goterm;
    }
    else
    {
      my @gotermArray;
      push @gotermArray, $goterm;
      $uniref2goterm{$uniref} = [ @gotermArray ];
    }
    if ($count % 1000000 == 0)
    {
      print STDERR "Read GOA: $count/132425703\r";
    }
  }
  close GOA;

#  print "========== GOA =========\n";
#  foreach my $uniref (keys %uniref2goterm)
#  {
#    print "$uniref: ";
#    foreach my $i (0 .. $#{ $uniref2goterm{$uniref} })
#    {
#      print  " $i = $uniref2goterm{$uniref}[$i]";
#    }
#    print "\n";
#  }

  # Create gene name, gene ontology term, and evalue.
  foreach my $locusid (keys %locus2protein)
  {
    my $protein = $locus2protein{$locusid};
    if (exists $protein2uniref{$protein})
    {
      my @unirefs = keys %{ $protein2uniref{$protein} };
      foreach my $uniref (@unirefs)
      {
        my $evalue = $protein2uniref{$protein}{$uniref};
        foreach my $i (0 .. $#{ $uniref2goterm{$uniref} })
        {
          print $outfile "$locusid\t$protein\t$uniref\t$uniref2goterm{$uniref}[$i]\t$evalue\n";
        }
      }
    }
  }
}
elsif ($cmd eq "go2ngene")
{
  # Parse gene2go.
  my @genes;
  open GENE2GO, $params{gene2go} or die "cannot open < $params{gene2go} $!";
  while (<GENE2GO>)
  {
    chomp;
    my $rec = {};
    ($rec->{gene}, $rec->{protein}, $rec->{uniref}, $rec->{goid},$rec->{evalue}) = split /\t/; 
    push @genes, $rec;
  }
  close GENE2GO;

  # Parse obo.
  my $line;
  my %goterms;
  open OBO, $params{obo} or die "cannot open < $params{obo} $!";
  while ($line = <OBO>)
  {
    chomp $line;
    if ($line =~ /^\[Term\]/)
    {
      $line = <OBO>; chomp $line;
      $line =~ /^id\:\s+([\w\:]+)$/;
      my $goid = $1;
      $line = <OBO>; chomp $line;
      $line =~ /^name\:\s+(.+)$/;
      my $goname = $1;
      $goterms{$goid} = $goname;
    }
  }
  close OBO;

  # Count genes for gene ontology ID.
  my %gocat;
  foreach my $i (0 .. $#genes)
  {
    my $g = $genes[$i];
    unless (exists $gocat{$g->{goid}})
    {
      my $rec = {};
      $rec->{ngene} = 0;
      $rec->{desc} = $goterms{$g->{goid}};
      $gocat{$g->{goid}} = $rec;
    }
    my $rec = $gocat{$g->{goid}}; 
    $rec->{ngene}++;
  }

  foreach my $key (keys %gocat)
  {
    print $outfile "$key\t$gocat{$key}{ngene}\t$gocat{$key}{desc}\n";
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

geneontology.pl - Functional category analyses

=head1 VERSION

geneontology.pl 1.0

=head1 SYNOPSIS

perl geneontology.pl gene2go -gff f.gff -blast f.blast -goa gene_association.goa_uniprot_noiea

perl geneontology.pl go2ngene -gene2go smutans.gene2go -obo gene_ontology.1_2.obo

=head1 DESCRIPTION

geneontology.pl produces files that are needed to functional category analyses.
One file is created using command gene2go: gene, gene ontology term, and real
number. Files can be created using command go2ngene: gene ontology term, number
of genes, and description of the term.

=head1 OPTIONS

  command: 

  gene2go - A file with tab-delimited elements is created using bacterial genome
  gff file, a tabular-form of BLAST output file, and goa file. The gff file is
  used to find protein IDs of gene names: e.g.,
  locus_tag=SMU.01;...;protein_id=NP720488.1. BLAST output file is used to find
  uniprotein ID for protein id: e.g., NP_720488.1 ... UniRef90_Q8DWN9. The goa
  file is used to find the gene ontology term of a uniprotein name: e.g.,
  GO:0005515 ... UniProtKB:Q9XEK5.

  go2ngene - A file with tab-dellimited elements is created using an output file
  from command gene2go, and gene ontology obo file: gene ontology term, number
  of genes, and description of the term. For each gene ontology term I count
  genes using an output file from command gene2go, and find the description of
  the term in an obo file.

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-gff> <file>

A bacterial genome gff annotation file.

=item B<-blast> <file>

A BLAST output file. The file must be created blastp option -outfmt 6.

=item B<-goa> <file>

A goa file: e.g.,

http://www.geneontology.org/gene-associations/submission/gene_association.goa_uniprot.gz

=item B<-obo> <file>

An obo file: e.g.,

http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo

=item B<-out> <file>

An output file.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make geneontology.pl better.

=head1 COPYRIGHT

Copyright (C) 2011 Sang Chul Choi

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
