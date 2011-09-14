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
            'reference=s',
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
  while (readdir $dh) 
  {
    next if $_ eq "." or $_ eq "..";
    my $isFirstSequence = 1;
    my $f = "$indir/$_";
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
          print $outfile "s\t$id\t$start\t$lenAlignment\t$strand\t$len\t$seq\n";
        }
        my $lineRightDirectRemoved = substr $line, 1;
        my @e = split /\s+/, $lineRightDirectRemoved;
        $id = $e[0];
        $len = $e[1];
        $start = $e[2];
        $end = $e[3];
        if ($start < $end)
        {
          $strand = '+';
          $start = $start - 1;
        }
        else
        {
          $strand = '-';
          $start = $len - $end;
        }
        $seq = "";
        if ($isFirstSequence == 1)
        {
          $isFirstSequence = 0;
          print $outfile "a\n";
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
    print $outfile "s\t$id\t$start\t$lenAlignment\t$strand\t$len\t$seq\n";
    print $outfile "\n";
    close ALN;
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

=head1 DESCRIPTION

extract-msa.pl parses a BLAST result file to create a list of FASTA files.

=head1 OPTIONS

  command: 

  blast6 - Read a BLAST result to create FASTA-format files.
  muscle - Convert FASTA-format alignments to MAF-format ones.

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
