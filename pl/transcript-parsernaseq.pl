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
my $VERSION = 'transcript-parsernaseq.pl 1.0';

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
            'parsernaseq=s',
            'wiggle=s',
            'feature=s',
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
my $parsernaseq;

if (exists $params{in})
{
  $in = $params{in};
}

if (exists $params{parsernaseq})
{
  $parsernaseq = $params{parsernaseq};
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

if ($cmd eq "pileup")
{
  unless (exists $params{wiggle})
  {
    &printError("Command $cmd needs -wiggle options");
  }
}
elsif ($cmd eq "gene")
{
  unless (exists $params{feature})
  {
    &printError("Command $cmd needs -feature options");
  }
}
elsif ($cmd eq "bed" or $cmd eq "gff")
{
  unless (exists $params{parsernaseq})
  {
    &printError("Command $cmd needs -parsernaseq options");
  }
}
elsif ($cmd eq "operon")
{
  unless (exists $params{parsernaseq} and exists $params{feature})
  {
    &printError("Command $cmd needs -parsernaseq and -feature options");
  }
}

###############################################################################
# DATA PROCESSING
###############################################################################
if ($cmd eq "pileup")
{
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

  for (my $i = 0; $i <= $#wx; $i++)
  {
    my $pos = $i + 1;
    print $outfile "$pos\t$wx[$i]\n";
  }
}
elsif ($cmd eq "gene")
{
  ##############################################################
  # Read the gene file.
  my @genes = rnaseqPosParse ($feature); 

  ##############################################################
  # Find genes and their coverages.
  for (my $i = 0; $i <= $#genes; $i++)
  {
    my $g = $genes[$i];
    my $l = $g->{end} - $g->{start} + 1;
    my $spaces = "    "; # ParseRNAseq needs spaces not tabs.
    print $outfile "$g->{name}$spaces$g->{strand}$spaces$g->{start}$spaces$g->{end}$spaces$l\n";
  }
}
elsif ($cmd eq "bed")
{
  open PARSERNASEQ, $parsernaseq or die "cannot open < $parsernaseq $!";
  # have a state
  my @m;
  my $i = 0;
  while (<PARSERNASEQ>)
  {
    $i++;
    my @e = split /\s+/;
    my $rec = {};
    $rec->{coverage} = $e[1];
    my $state = $e[2];
    if ($state =~ /C(\d+)/)
    {
      $rec->{state} = $1;
    }
    else
    {
      $rec->{state} = 0;
    }
    push @m, $rec;
    if ($i % 100000 == 0)
    {
      print STDERR "Reading $i\r";
    }
  }
  close PARSERNASEQ;

  # Find 
  my @bd;
  my $prev = "";
  my %operon;
  for (my $i = 0; $i <= $#m; $i++)
  {
    my $rec = $m[$i];
    my $curr = $rec->{state};
    if ($curr eq $prev)
    {
    }
    else
    {
      unless ($i == 0)
      {
        $operon{end} = $i + 1;
        push @bd, { %operon };
      }
      $operon{start} = $i + 1;
      $operon{state} = $curr;
    }
    $prev = $curr;
  }

  # Print
  for (my $i = 0; $i <= $#bd; $i++)
  {
    my $rec = $bd[$i];
    print $outfile "chr1\t$rec->{start}\t$rec->{end}\t$rec->{state}\t$rec->{state}\t+\t$rec->{start}\t$rec->{end}\t255,0,0\n";
  }
}
elsif ($cmd eq "gff")
{
  open PARSERNASEQ, $parsernaseq or die "cannot open < $parsernaseq $!";
  while (<PARSERNASEQ>)
  {
    chomp;
    my @e = split /\t/;
    my $start = $e[3];
    my $end   = $e[4];
    my $v1    = $e[5];
    my $v2    = $e[8];
    $v2 =~ /Name=(C\d+)\_/; 
    my $v3 = $1;
    print $outfile "chr1\t$start\t$end\t$v3\t$v1\t+\n";
  }
  close PARSERNASEQ;
}
elsif ($cmd eq "operon")
{
  ##############################################################
  # Read the gene file.
  my @genes = rnaseqPosParse ($feature); 

  #################################################################
  # Read ParseRNAseq in GFF format.
  my @m;
  open PARSERNASEQ, $parsernaseq or die "cannot open < $parsernaseq $!";
  while (<PARSERNASEQ>)
  {
    chomp;
    my @e = split /\t/;
    my $rec = {};
    $rec->{start} = $e[3];
    $rec->{end}   = $e[4];
    $rec->{score} = $e[5];
    $e[8] =~ /Name=(.+)/; 
    $rec->{name}  = $1;
    $rec->{strand} = '+';
    push @m, $rec;
  }
  close PARSERNASEQ;

  #################################################################
  # Group genes within transcrips.
  for (my $i = 0; $i <= $#m; $i++)
  {
    my $t = $m[$i];
    my @genesT;
    for (my $j = 0; $j <= $#genes; $j++)
    {
      my $g = $genes[$j];
      if (($t->{start} <= $g->{start} and $g->{start} <= $t->{end})
          or ($t->{start} <= $g->{end} and $g->{end} <= $t->{end})
          or ($g->{start} <= $t->{start} and $t->{end} <= $g->{end}))
      {
        my $rec = {};
        $rec->{name}   = $g->{name};
        $rec->{start}  = $g->{start};
        $rec->{end}    = $g->{end};
        $rec->{strand} = $g->{strand};
        push @genesT, $rec;
      }
    }
    $t->{gene} = \@genesT;
  }

  for (my $i = 0; $i <= $#m; $i++)
  {
    my $t = $m[$i];
    print "$t->{name}: ";
    my @genesT;
    if (exists $t->{gene})
    {
      @genesT = @{$t->{gene}};
      for (my $j = 0; $j <= $#genesT; $j++)
      {
        my $g = $genesT[$j];
        print ":$g->{name}";
      } 
    }
    print "\n";
  }
  print "\n++++++++++++++++++++++++++++++++\n";

  my @m2;
  for (my $i = 0; $i <= $#m; $i++)
  {
    my $t = $m[$i];
    my $recT = {};
    $recT->{name}  = $t->{name};
    $recT->{start} = $t->{start};
    $recT->{end}   = $t->{end};
    $recT->{score} = $t->{score};
    $recT->{strand} = $t->{strand};
    push @m2, $recT;

    my @genesT = @{$t->{gene}};

    my @genesT2;
    my $g = $genesT[0];
    my $prevStrand = $g->{strand};
    for (my $j = 0; $j <= $#genesT; $j++)
    {
      $g = $genesT[$j];

      my $rec = {};
      $rec->{name}   = $g->{name};
      $rec->{start}  = $g->{start};
      $rec->{end}    = $g->{end};
      $rec->{strand} = $g->{strand};
      my $currStrand = $g->{strand};
      if ($j == 0 and $j < $#genesT)
      {
        push @genesT2, $rec;
        $prevStrand = $g->{strand};
      }
      elsif ($j == 0 and $j == $#genesT)
      {
        push @genesT2, $rec;
        $recT = $m2[$#m2];
        $recT->{gene} = \@genesT2;
      }
      elsif ($j == $#genesT)
      {
        if ($currStrand eq $prevStrand)
        {
          push @genesT2, $rec;
          $recT = $m2[$#m2];
          $recT->{gene} = \@genesT2;
        }
        else
        {
          $recT = $m2[$#m2];
          my @genesT2copy = @genesT2;
          $recT->{gene} = \@genesT2copy;

          $recT = {};
          $recT->{name}  = $t->{name};
          $recT->{start} = $t->{start};
          $recT->{end}   = $t->{end};
          $recT->{score} = $t->{score};
          $recT->{strand} = $t->{strand};
          @genesT2 = ();
          push @genesT2, $rec;
          my @genesT2copy2 = @genesT2;
          $recT->{gene}  = \@genesT2copy2;
          push @m2, $recT;
        }
      }
      else
      {
        if ($currStrand eq $prevStrand)
        {
          push @genesT2, $rec;
        }
        else
        {
          $recT = $m2[$#m2];
          my @genesT2copy = @genesT2;
          $recT->{gene} = \@genesT2copy;

          @genesT2 = ();
          push @genesT2, $rec;
          $recT = {};
          $recT->{name}  = $t->{name};
          $recT->{start} = $t->{start};
          $recT->{end}   = $t->{end};
          $recT->{score} = $t->{score};
          $recT->{strand} = $t->{strand};
          push @m2, $recT;
        }
      }
      $prevStrand = $currStrand;
    }
  }

  for (my $i = 0; $i <= $#m2; $i++)
  {
    my $t = $m2[$i];

    if (exists $t->{gene})
    {
      my @genesT = @{$t->{gene}};
      my $g = $genesT[0];
      $t->{start} = $g->{start};
      $g = $genesT[$#genesT];
      $t->{end} = $g->{end};
      $t->{strand} = $g->{strand};
    }
    print "\n";
  }

  for (my $i = 0; $i <= $#m2; $i++)
  {
    my $t = $m2[$i];
    print $outfile "chr1\t$t->{start}\t$t->{end}\t$t->{name}\t$t->{score}\t$t->{strand}\n";
    # print $outfile "chr1\t$t->{start}\t$t->{end}\t$t->{name}\t$t->{score}\t$t->{strand}\n";
  }
}

if (exists $params{out})
{
  close $outfile;
}
__END__
=head1 NAME

transcript-parsernaseq - input and output files for ParseRNASeq

=head1 VERSION

transcript-parsernaseq 1.0

=head1 SYNOPSIS

perl pl/transcript-parsernaseq.pl pileup -wiggle FASTQ01.wig -out FASTQ01.parsernaseq.pileup

perl pl/transcript-parsernaseq.pl gene -feature feature-genome.out-geneonly -out out.parsernaseq

perl pl/transcript-parsernaseq.pl bed -parsernaseq out.parsernaseq -out bed.file

perl pl/transcript-parsernaseq.pl gff -parsernaseq out.parsernaseq -out bed.file

perl pl/transcript-parsernaseq.pl operon -feature feature-genome.out-geneonly -parsernaseq out.parsernaseq -out operon.file

=head1 DESCRIPTION

transcript-parsernaseq will help you to pre- and post-process files of
ParseRNASeq.

Command:
  pileup - Using FASTQ01.wig file that is created by samtools mpileup to create
  a pileup file for ParaseRNASeq.

  gene - Using feature-genome.out-geneonly file that is created by menu
  feature-genome I create a gene annotation file for ParseRNASeq.

  bed - Using ParseRNAseq output to create a bed file.

  gff - Using ParseRNAseq gff output to create a bed file.

  operon - Using ParseRNAseq gff output to create an operon file.

=head1 OPTIONS

=over 8

=item B<-feature> <file>

Features of a genome. 

=item B<-wiggle> <file>

The wiggle file of coverage values.

=item B<-parsernaseq> <file>

The output file from ParseRNAseq.

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
