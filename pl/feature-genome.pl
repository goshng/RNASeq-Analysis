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
my $VERSION = 'feature-genome.pl 1.0';

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
            'in=s',
            'chromosome=s',
            'bed=s',
            'feature=s',
            'intergenicregion',
            'intergenicregiononly',
            'startcodon',
            'geneonly',
            'withproduct',
            'windowsize=i',
            'genomelength=i',
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
my $feature;
my $chromosome = "chr1";
my $minLenNoncoding = 50;
my $bed;
my $windowsize = 200;
my $genomelength;

if (exists $params{bed}) {
  $bed = $params{bed};
}

if (exists $params{genomelength}) {
  $genomelength = $params{genomelength};
}

if (exists $params{windowsize}) {
  $windowsize = $params{windowsize};
}

if (exists $params{chromosome}) {
  $chromosome = $params{chromosome};
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

if ($cmd eq "extract")
{
  unless (exists $params{bed})
  {
    &printError("command extract needs option -bed");
  }
}

################################################################################
## DATA PROCESSING
################################################################################
sub convert_gff_ingene($$$); 

if ($cmd eq "ptt")
{
  my @igr;
  my $prevEnd = 0;
  my $line;
  while ($line = <$infile>)
  {
    next unless $line =~ /^(\d+)\.\.(\d+)\s+/;
    my @e = split /\t/, $line;
    $e[0] =~ /^(\d+)\.\.(\d+)\s+/;
    my $start  = $1 - 1; # PTT is 1-base, BED is 0-base.
    my $end    = $2;
    my $strand = $e[1];
    my $name   = $e[5];

    my $lengthNoncoding = $start - $prevEnd;
    if ($lengthNoncoding > $minLenNoncoding)
    {
      my $startNoncoding = $prevEnd; # BED is 0-base
      my $startNoncodingOnebase = $prevEnd + 1;
      my $noncoding = {};
      $noncoding->{name}   = "IGR$startNoncodingOnebase-$start";
      $noncoding->{start}  = $startNoncoding;
      $noncoding->{end}    = $start;
      $noncoding->{strand} = '+';
      push @igr, $noncoding;
      if (exists $params{intergenicregion}
          or exists $params{intergenicregiononly})
      {
        print $outfile "$chromosome\t";
        print $outfile "$noncoding->{start}\t";
        print $outfile "$noncoding->{end}\t";
        print $outfile "$noncoding->{name}\t0\t";
        print $outfile "$noncoding->{strand}\n";
      }
    }
    $prevEnd = $end;
    unless (exists $params{intergenicregiononly})
    {
      print $outfile "$chromosome\t$start\t$end\t$name\t0\t$strand\n";
    }
  }
}
elsif ($cmd eq "ptt2")
{
  # Read all of the PTT file.
  my @ptt = rnaseqPttParse ($in);
  $genomelength = rnaseqPttGetGenomeLength ($in); 

  # Then, create a BED file depending on options.
  if (exists $params{startcodon})
  {
    for (my $i = 0; $i <= $#ptt; $i++)
    {
      my $g = $ptt[$i];
      $g->{Location} =~ /(\d+)\.\.(\d+)/;
      my $start = $1 - 1; # PTT is 1-base, BED is 0-base.
      my $end = $2;
      if ($g->{Strand} eq '+')
      {
        $end    = $start + $windowsize + 1;
        $start  -= $windowsize;
      }
      else
      {
        $start = $end - $windowsize - 1;
        $end   += $windowsize;
      }
      if ($start < 0)
      {
        $start = 0;
      }
      if ($end > $genomelength)
      {
        $end = $genomelength;
      }

      if (exists $params{withproduct})
      {
        print $outfile "$chromosome\t$start\t$end\t$g->{Synonym}\t0\t$g->{Strand}\t";
        print $outfile "$g->{Product}\n";
      }
      else
      {
        print $outfile "$chromosome\t$start\t$end\t$g->{Synonym}\t0\t$g->{Strand}\n";
      }
    }
  }
  elsif (exists $params{geneonly})
  {
    for (my $i = 0; $i <= $#ptt; $i++)
    {
      my $g = $ptt[$i];
      $g->{Location} =~ /(\d+)\.\.(\d+)/;
      my $start = $1 - 1; # PTT is 1-base, BED is 0-base.
      my $end = $2;

      if (exists $params{withproduct})
      {
        print $outfile "$chromosome\t$start\t$end\t$g->{Synonym}\t0\t$g->{Strand}\t";
        print $outfile "$g->{Product}\n";
      }
      else
      {
        print $outfile "$chromosome\t$start\t$end\t$g->{Synonym}\t0\t$g->{Strand}\n";
      }
    }
  }
}
elsif ($cmd eq "gff")
{
  convert_gff_ingene($in, $feature, $out); 
}
elsif ($cmd eq "extract")
{
  # I might want to use twoBitTwoFa.
  my @bedfeature = rnaseqBedParse ($bed);
  my $seq = maFastaParse ($in);
  for (my $i = 0; $i <= $#bedfeature; $i++)
  {
    my $b = $bedfeature[$i];
    my $offset = $b->{start};
    my $l = $b->{end} - $b->{start};
    my $subseq = substr $seq, $offset, $l;
    if (exists $b->{strand})
    {
      if ($b->{strand} eq '-')
      {
        $subseq =~ tr/AGCTUagctu/TCGAAtcgaa/;
        $subseq = reverse $subseq;
      }
    }
    print $outfile ">$b->{name}\n";
    print $outfile "$subseq\n";
  }
}
elsif ($cmd eq "parsernaseq")
{
  my $i = 0;
  while (<$infile>)
  {
    $i++;
    my $istr = sprintf ("%04d", $i);
    my @e = split /\t/;
    print $outfile "$chromosome\t";
    print $outfile "$e[3]\t$e[4]\ttx$istr\t0\t+\n";
  }
}
elsif ($cmd eq "rnaz")
{
  while (<$infile>)
  {
    chomp;
    print $outfile "$_+\n";
  }
}
elsif ($cmd eq "keg")
{
  my %c1;
  my $l;
  my $count = 0;
  my $count2 = 0;
  my $level1 = "";
  my $level2 = "";
  my $level3 = "";
  my $level4 = "";
  my @genes = ();
  my %categoryDescription;
  my %categoryCountDescription;
  while ($l = <$infile>)
  {
    chomp $l;
    unless ($l =~ /^D/)
    {
      if (scalar (@genes) > 0)
      {
        $count2++;
        foreach my $g (@genes) 
        {
          $g =~ s/\_/\./g;
          print $outfile "$level1\t$level2\t$level3\t$g\n"
        }
        my $rec = {};
        $categoryCountDescription{$level3} = $rec;
        $rec->{description} = $categoryDescription{$level3};
        $rec->{count} = scalar (@genes);
        #print "$count2 (1) $level1 : (2) $level2 : (3) $level3\t";
        #print (join(", ", @genes));
        #print "\n";
        @genes = ();
      }
    }

    if ($l =~ /^A\<b\>(.+)\<\/b\>/)
    {
      $level1 = $1;
      # print $level1,"\n";
    }
    elsif ($l =~ /^B\s+\<b\>(.+)\<\/b\>/)
    {
      $level2 = $1;
      # print "\t",$level2,"\n";
    }
    # elsif ($l =~ /^C\s+(\d+)\s+(.+)\s+\[PATH:smu\d+\]/)
    elsif ($l =~ /^C\s+(\d+)\s+(.+)/)
    {
      $level3 = "KEGG:$1";
      $count++;
      # print "\t\t",$level3,"\n";
      $categoryDescription{$level3} = $2;
    }
    elsif ($l =~ /^D\s+(\w+)/)
    {
      $level4 = $1;
      # print "\t\t\t",$level4,"\n";
      push @genes, $level4;
    }
#        00010 Glycolysis / Gluconeogenesis [PATH:smu00010]

#      $rec = {};
#      $c1{$1} = $rec;
#      for $field ( split ) {
#        ($key, $value) = split /=/, $field;
#        $rec->{$key} = $value;
#      }

  }

  # print "Count: $count\n";

  # Print $categoryDescription{$1} = $2;
  foreach my $cat (keys %categoryCountDescription)
  {
    print $outfile "$cat\t$categoryCountDescription{$cat}{count}";
    print $outfile "\t$categoryCountDescription{$cat}{description}\n";
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

sub convert_gff_ingene($$$) {
  my ($gffFilename, $f, $out) = @_;
  open OUT, ">$out" or die "$out could not be opened"; 

  my $numGene = 0;
  my $gffLine;
  my $blockid_prev = 0;
  open GFF, "$gffFilename" or die $!;
  while ($gffLine = <GFF>)
  {
    if ($feature eq "gene")
    {
      if ($gffLine =~ /RefSeq\s+gene\s+(\d+)\s+(\d+).+([+-]).+locus_tag=([\w\.]+);/)
      {
        print OUT "$chromosome\t$1\t$2\t$3\t$4\n";
      }
    }
    elsif ($feature eq "genestart")
    {
      if ($gffLine =~ /RefSeq\s+gene\s+(\d+)\s+(\d+).+([+-]).+locus_tag=([\w\.]+);/)
      {
        print OUT "$chromosome\t$1\t$3\n";
      }
    }

  }
  close GFF;
  close OUT;
}
__END__
=head1 NAME

feature-genome.pl - Parse a gff file to make ingene file

=head1 VERSION

feature-genome.pl 1.0

=head1 SYNOPSIS

perl feature-genome.pl ptt -in file.ptt

perl feature-genome.pl ptt -in file.ptt -intergenicregion

perl feature-genome.pl ptt -in file.ptt -startcodon

perl feature-genome.pl ptt2 -startcodon -windowsize 100 -in file.ptt 

perl pl/feature-genome.pl keg -in data/smu00001.keg

=head1 DESCRIPTION

feature-genome.pl parses a genome annotation file such as a gff file. A BED output
file is created with tab-delimited columns: chrmosome, start, end, name, strand.
The minimum length of intergenic regions is 30 base pairs long. The BED file can
be used to create a FASTA-format file with DNA sequences that are extracted from
a reference genome. Note that BED files have 0-base start positions. The length
of a feature should be obtained by subtracting the start position from the end
position.

=head1 OPTIONS

  command: 

  ptt - Reads a ptt file to create a BED file.

  gff - Reads a gff file.

  extract - Extract DNA sequences using a BED file.
  I extract a reverse complementary sequences for the regions with minus strand.

  ptt2 - Use a ptt file to create a BED file. This might replace the command
  ptt. With option -startcodon regions centered around the start codon of genes
  are extracted. Use option -windowsize to change the size of regions. Default
  of -windowsize is 200: The length of a block is 401.

  keg - use a keg file to create a gene2cat file. Use a sample file, data/smu00001.keg 
  The contains three levels of functional categories. The third level contains
  genes. Note that gene names specification may be different from that of
  GenBank. The output file cotains two two parts: the first part one is of
  four-column, and the second part is two-column. The first, second, and third
  columns of the first part are the three levels of KEGG pathways. The fourth
  column contains gene names. The second part of the file contains descriptions
  of the third level, or the third column of the first part. These two parts
  could be used in functional category association.

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

=item B<-out> <file>

An output file.

=item B<-intergenicregion>

Not only genes but also intergenic regions are also added to the output file.

=item B<-bed>

A BED file.

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make feature-genome.pl better.

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
