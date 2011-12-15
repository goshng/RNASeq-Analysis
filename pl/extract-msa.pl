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
$| = 1; # Do not buffer output
my $VERSION = 'extract-msa.pl 1.0';

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
            'muscle=s',
            'reference=s',
            'genomeLength=i',
            'in=s',
            'indir=s',
            'out=s',
            'outdir=s',
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
my $outdir;
my $indir;
my $reference;
my $muscle;
my $genomeLength;

if (exists $params{genomeLength}) {
  $genomeLength = $params{genomeLength};
}

if (exists $params{muscle}) {
  $muscle = $params{muscle};
}

if (exists $params{reference}) {
  $reference = $params{reference};
}

if (exists $params{indir}) {
  $indir = $params{indir};
}

if (exists $params{outdir}) {
  $outdir = $params{outdir};
}

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

if ($cmd eq "blast6")
{
  unless (exists $params{outdir}
          and exists $params{reference})
  {
    &printError("command blast6 needs options -outdir and -reference");
  }
}
elsif ($cmd eq "muscle")
{
  unless (exists $params{indir})
  {
    &printError("command muscle needs options -indir and -in");
  }
}
elsif ($cmd eq "clust")
{
  unless (exists $params{in}
          and exists $params{reference}
          and exists $params{muscle}
          and exists $params{outdir})
  {
    &printError("command muscle needs options -indir and -in");
  }
}

################################################################################
## DATA PROCESSING
################################################################################
if ($cmd eq "blast6")
{
  my $qseqid = "first";
  my @sequences;
  my @sequenceIDs;
  my $alnNamePrev = "first";
  my $line;
  while ($line = <$infile>)
  {
    chomp $line;
    my @e = split /\t/, $line;
    my $h = {};
    $h->{qseqid}   = $e[0]; 
    $h->{sseqid}   = $e[1];
    $h->{pident}   = $e[2]; 
    $h->{length}   = $e[3]; 
    $h->{mismatch} = $e[4]; 
    $h->{gapopen}  = $e[5]; 
    $h->{qstart}   = $e[6]; 
    $h->{qend}     = $e[7]; 
    $h->{sstart}   = $e[8]; 
    $h->{send}     = $e[9]; 
    $h->{evalue}   = $e[10]; 
    $h->{bitscore} = $e[11]; 
    $h->{qlen}     = $e[12]; 
    $h->{qseq}     = $e[13]; 
    $h->{slen}     = $e[14]; 
    $h->{sseq}     = $e[15]; 
    my $id = "$h->{sseqid}-$h->{sstart} $h->{slen} $h->{sstart} $h->{send}";
    if ($alnNamePrev eq $h->{qseqid})
    {
      # Requirements for orthologs.
      # 1. It must not be from the reference genome.
      # 2. At least 70% of the total length of the query sequence.
      if ($h->{pident} == 100 and $h->{qlen} == $h->{length})
      {
        # Do not include this duplcate sequence.
      }
      else
      {
        unless ($h->{sseqid} =~ /$reference/)
        {
          if ($h->{length} > 0.7 * $h->{qlen})
          {
            push @sequenceIDs, $id;
            $h->{sseq} =~ s/-//g;
            push @sequences, $h->{sseq};
          }
        }
      }
    }
    else
    {
      unless ($qseqid eq "first")
      {
        # Close the previous file.
        if ($#sequences > 0)
        {
          open FASTA, ">", $qseqid or die "cannot open > $qseqid $!";
          for (my $i = 0; $i <= $#sequences; $i++)
          {
            print FASTA ">$sequenceIDs[$i]\n";
            print FASTA "$sequences[$i]\n";
          }
          close FASTA;
        }
      }
      # Create a new file.
      $qseqid = "$outdir/$h->{qseqid}.fa";
      @sequenceIDs = ();
      @sequences = ();
      push @sequenceIDs, $id;
      push @sequences, $h->{sseq};
      unless ($h->{sseqid} =~ /$reference/)
      {
        die "The first sequence is not reference";
      }
    }
    $alnNamePrev = $h->{qseqid};
  }
  # Last entry.
  if ($#sequences > 0)
  {
    open FASTA, ">", $qseqid or die "cannot open > $qseqid $!";
    for (my $i = 0; $i <= $#sequences; $i++)
    {
      print FASTA ">$sequenceIDs[$i]\n";
      print FASTA "$sequences[$i]\n";
    }
    close FASTA;
  }
}
elsif ($cmd eq "clust")
{
  #######################################
  # Find lines with the first alignments.
  # @linenumbersFirstAlignment contains start line numbers of queries, and the
  # number of lines in the BLAST file.
  my @linenumbersFirstAlignment;
  my $alnNamePrev = "first";
  my $line;
  my $n = 0;
  while ($line = <$infile>)
  {
    $n++;
    my @e = split /\t/, $line;
    my $h = {};
    $h->{qseqid}   = $e[0]; 
    unless ($alnNamePrev eq $h->{qseqid})
    {
      push @linenumbersFirstAlignment, $n;
    }
    $alnNamePrev = $h->{qseqid};
  }
  $n++;
  push @linenumbersFirstAlignment, $n;
  close $infile;

  ###############################
  # Start of the clust operation.
  open ($infile, "<", $in) or die "cannot open < $in: $!";
  for (my $i = 0; $i < $#linenumbersFirstAlignment; $i++)
  {
    my $numberFailedToBeAdded = 0;
    my @cluster;
    my $lineStart = $linenumbersFirstAlignment[$i];
    my $lineEnd = $linenumbersFirstAlignment[$i+1];
    my $fastaFilename; 
    for (my $j = 0; $j < $lineEnd - $lineStart; $j++)
    {
      $line = <$infile>;
      chomp $line;
      my @e = split /\t/, $line;
      my $h = {};
      $h->{qseqid}   = $e[0]; 
      $h->{sseqid}   = $e[1];
      $h->{pident}   = $e[2]; 
      $h->{length}   = $e[3]; 
      $h->{mismatch} = $e[4]; 
      $h->{gapopen}  = $e[5]; 
      $h->{qstart}   = $e[6]; 
      $h->{qend}     = $e[7]; 
      $h->{sstart}   = $e[8]; 
      $h->{send}     = $e[9]; 
      $h->{evalue}   = $e[10]; 
      $h->{bitscore} = $e[11]; 
      $h->{qlen}     = $e[12]; 
      $h->{qseq}     = $e[13]; 
      $h->{slen}     = $e[14]; 
      $h->{sseq}     = $e[15]; 
      if ($j == 0)
      {
        $fastaFilename = "$outdir/$h->{qseqid}.fa";
      }

      # 1. Check the first alignment.
      if ($j == 0)
      {
        unless ($h->{pident} == 100 and $h->{length} == $h->{qlen})
        {
          die "Each query sequence must be aligned the reference genome";
        }
      }
      else
      {
        # Filter sequences out: the subject sequence from the reference, 
        # genomes that are already added. I also skip to try to align sequences
        # if the number of aligning attempts exceeds 10.
        my $hc;
        if ($h->{sseqid} =~ /$reference/)
        {
          next;
        }
        my $isDuplicate = 0;
        for (my $ci = 0; $ci <= $#cluster; $ci++)
        {
          $hc = $cluster[$ci];
          if ($hc->{sseqid} eq $h->{sseqid})
          {
            $isDuplicate = 1;
            last;
          }
        }
        if ($isDuplicate == 1)
        {
          next;
        }
        if (scalar @cluster >= 1)
        {
          $hc = $cluster[0];
          $h->{sseq} =~ s/-//g;
          if (length ($h->{sseq}) < 0.5 * $hc->{length})
          {
            next;
          }
        }
        if ($numberFailedToBeAdded > 10)
        {
          next;
        }

        # Create a FASTA-format file with two sequences.
        if (scalar @cluster == 0)
        {
          rnaseqFastaCreate ($fastaFilename, "chr1", $h->{qseq});
          rnaseqFastaAdd ($fastaFilename, $h->{sseqid}, $h->{sseq});
          push @cluster, $h;
        }
        elsif (scalar @cluster >= 1)
        {
          $hc = $cluster[0];
          # Create the final FASTA file for the current query.
          for (my $ci = 0; $ci <= $#cluster; $ci++)
          {
            $hc = $cluster[$ci];
            if ($ci == 0)
            {
              rnaseqFastaCreate ($fastaFilename, "chr1", $hc->{qseq});
              rnaseqFastaAdd ($fastaFilename, $hc->{sseqid}, $hc->{sseq});
            }
            else
            {
              rnaseqFastaAdd ($fastaFilename, $hc->{sseqid}, $hc->{sseq});
            }
          }
          rnaseqFastaAdd ($fastaFilename, $h->{sseqid}, $h->{sseq});

          # Align the sequences using MUSCLE.
          my $muscleCommand = "$muscle -maxiters 1 -diags -quiet -in $fastaFilename -out $fastaFilename.muscle";
          system ($muscleCommand);

          # Keep the current sequence if the alignment is good.
          my $r = rnaseqFastaManyGap ("$fastaFilename.muscle");
          if ($r == 0)
          {
            push @cluster, $h;
            $numberFailedToBeAdded = 0; 
          }
          else
          {
            $numberFailedToBeAdded++; 
          }
        }
      }
    }

    unlink "$fastaFilename";
    unlink "$fastaFilename.muscle";

    # Create the final FASTA file for the current query.
    if ($#cluster >= 1)
    {
      my $h = $cluster[0];
      for (my $j = 0; $j <= $#cluster; $j++)
      {
        $h = $cluster[$j];
        if ($j == 0)
        {
          rnaseqFastaCreate ($fastaFilename, "chr1 $h->{qlen} $h->{qstart} $h->{qend}", 
                             $h->{qseq});
          rnaseqFastaAdd ($fastaFilename, 
                          "$h->{sseqid} $h->{slen} $h->{sstart} $h->{send}", 
                          $h->{sseq});
        }
        else
        {
          rnaseqFastaAdd ($fastaFilename, 
                          "$h->{sseqid} $h->{slen} $h->{sstart} $h->{send}", 
                          $h->{sseq});
        }
      }
    }
    print STDERR "Processing $i/$#linenumbersFirstAlignment\r";
  }
}
elsif ($cmd eq "muscle")
{
  my $seq;
  my $line;
  my $id;
  my $len;
  my $start;
  my $end;
  my $strand;
  opendir (my $dh, $indir) or die "cannot open $indir";
  while (my $fileDh = readdir ($dh)) 
  {
    next if $fileDh eq "." or $fileDh eq "..";

    # The first sequence must be from the reference genome.
    my @alignedSequences;
    my $h = {};
    $h->{id} = "chr1";
    push @alignedSequences, $h;

    my $isFirstSequence = 1;
    $fileDh =~ /IGR(\d+)-(\d+)\.fa/;
    my $qStart = $1;
    my $qEnd = $2;
    my $f = "$indir/$fileDh";
    open ALN, $f or die "cannot open < $f $!";
    while ($line = <ALN>)
    {
      chomp $line;
      if ($line =~ /^>/)
      {
        if ($isFirstSequence == 0)
        {
          my $seqWithoutGap = $seq;
          $seqWithoutGap =~ s/-//g;
          my $lenAlignment = length $seqWithoutGap; 
          if ($id eq "chr1")
          {
            my $hFirst = $alignedSequences[0];
            $hFirst->{start} = $start + ($qStart - 1); # start (0-base)
                                                       # qStart (1-base)
            $hFirst->{lenAlignment} = $lenAlignment;
            $hFirst->{strand} = $strand;
            $hFirst->{len} = $genomeLength;
            $hFirst->{seq} = $seq;
          }
          else
          {
            my $h = {};
            $h->{id} = $id;
            $h->{start} = $start;
            $h->{lenAlignment} = $lenAlignment;
            $h->{strand} = $strand;
            $h->{len} = $len;
            $h->{seq} = $seq;
            push @alignedSequences, $h;
          }
        }
        my $lineRightDirectRemoved = substr $line, 1;
        my @e = split /\s+/, $lineRightDirectRemoved;
        $id    = $e[0];
        $len   = $e[1];
        $start = $e[2];
        $end   = $e[3];
        if ($start < $end)
        {
          $strand = '+';
          $start = $start - 1;
        }
        else
        {
          $strand = '-';
          # $start = $len - $end;
          $start = $end - 1; # The end is the start.
        }
        $seq = "";
        if ($isFirstSequence == 1)
        {
          $isFirstSequence = 0;
        }
      }
      else
      {
        $seq .= $line;
      }
    }
    # Print the last sequence.
    my $seqWithoutGap = $seq;
    $seqWithoutGap =~ s/-//g;
    my $lenAlignment = length $seqWithoutGap; 
    if ($id eq "chr1")
    {
      my $hFirst = $alignedSequences[0];
      $hFirst->{start} = $start + ($qStart - 1); # start (0-base)
                                                 # qStart (1-base)
      $hFirst->{lenAlignment} = $lenAlignment;
      $hFirst->{strand} = $strand;
      $hFirst->{len} = $genomeLength;
      $hFirst->{seq} = $seq;
    }
    else
    {
      my $h = {};
      $h->{id} = $id;
      $h->{start} = $start;
      $h->{lenAlignment} = $lenAlignment;
      $h->{strand} = $strand;
      $h->{len} = $len;
      $h->{seq} = $seq;
      push @alignedSequences, $h;
    }
    close ALN;

    print $outfile "a\n";
    for (my $i = 0; $i <= $#alignedSequences; $i++)
    {
      my $h = $alignedSequences[$i];
      if ($h->{id} =~ /gi\|\d+\|ref\|(\w+)\.\d+\|/)
      {
        $h->{id} = $1;
      }
      print $outfile "s\t$h->{id}\t$h->{start}\t$h->{lenAlignment}";
      print $outfile "\t$h->{strand}\t$h->{len}\t$h->{seq}\n";
    }
    print $outfile "\n";
  }
  closedir $dh;
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

extract-msa.pl - Multiple sequences are extracted from BLAST results.

=head1 VERSION

extract-msa.pl 1.0

=head1 SYNOPSIS

perl extract-msa.pl blast6 -in file.blast -outdir outdir

perl extract-msa.pl muscle -indir indir -in file.blast -outdir outdir

perl extract-msa.pl clust -in file.blast -muscle muscle -outdir outdir

=head1 DESCRIPTION

extract-msa.pl parses a BLAST result file to create a list of FASTA files.

=head1 OPTIONS

  command: 

  blast6 - Read a BLAST result to create FASTA-format files.
  muscle - Convert FASTA-format alignments to MAF-format ones.
  clust - Find similar sequences for each intergenic region.

clust: for each intergenic region do the following steps. 1. Check the first
alignment if the query sequence is aligned to a position uniquely, 2. Create a
FASTA-format file with two sequences from following alignments; the subject
sequence must be from a new species (and it must not be from the reference
genome), 3. Add a new subject sequence from another following alignment to the
FASTA-format file, 4. Align the sequences in the FASTA file, 5. Keep the last
sequence if it is aligned very well (or it does not have 40 contiguous gaps), 6.
Go back to step 3 and repeat the procedure of alignmnt quality until there are
no alignments left.

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-in> <file>

A input file name.

=item B<-indir> <directory>

An input directory.

=item B<-outdir> <directory>

An output directory.

=item B<-out> <file>

An output file.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make extract-msa.pl better.

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
