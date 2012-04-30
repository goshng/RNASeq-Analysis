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
require 'pl/sub-array.pl';
require 'pl/sub-pos.pl';
require 'pl/sub-fasta.pl';
require 'pl/sub-wiggle.pl';
require 'pl/sub-bed.pl';

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
            'operon=s',
            'wiggle=s',
            'feature=s',
            'end=s',
            'inpileup=s',
            'peakcutoff=i',
            'sizecutoff=i',
            'windowsize=i',
            'toplevel',
            'updowncutoff=i',
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
my $end;
my $operon;
my $peakcutoff = 0.3;
my $sizecutoff = 70;
my $updowncutoff = 200;
my $windowsize = 100;

if (exists $params{in})
{
  $in = $params{in};
}

if (exists $params{windowsize})
{
  $windowsize = $params{windowsize};
}

if (exists $params{updowncutoff})
{
  $updowncutoff = $params{updowncutoff};
}

if (exists $params{peakcutoff})
{
  $peakcutoff = $params{peakcutoff};
}

if (exists $params{sizecutoff})
{
  $sizecutoff = $params{sizecutoff};
}

if (exists $params{operon})
{
  $operon = $params{operon};
}

if (exists $params{end})
{
  $end = $params{end};
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
elsif ($cmd eq "bed" or $cmd eq "gff" or $cmd eq "combine")
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
elsif ($cmd eq "adjust" or $cmd eq "slope")
{
  unless (exists $params{operon} and exists $params{end})
  {
    &printError("Command $cmd needs -operon and -end options");
  }
}
elsif ($cmd eq "remove")
{
  unless (exists $params{parsernaseq} and 
          exists $params{inpileup} and 
          exists $params{feature} and
          exists $params{out})
  {
    &printError("Command $cmd needs -parsernaseq, -feature and -inpileup options");
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
    my $l = $g->{end} - $g->{start};
    $g->{start}++;
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
  # Find the highest level.
  my $topLevel = 1;
  open PARSERNASEQ, $parsernaseq or die "cannot open < $parsernaseq $!";
  while (<PARSERNASEQ>)
  {
    chomp;
    my @e = split /\t/;
    my $v2    = $e[8];
    $v2 =~ /Name=C(\d+)\_/; 
    my $v3 = $1;
    if ($topLevel < $v3)
    {
      $topLevel = $v3;
    }
  }
  close PARSERNASEQ;

  open PARSERNASEQ, $parsernaseq or die "cannot open < $parsernaseq $!";
  while (<PARSERNASEQ>)
  {
    chomp;
    my @e = split /\t/;
    my $start = $e[3];
    my $end   = $e[4];
    my $v1    = $e[5];
    my $v2    = $e[8];
    $v2 =~ /Name=C(\d+)\_/; 
    my $v3 = $1;
    if (exists $params{toplevel})
    {
      if ($v3 == $topLevel)
      {
        print $outfile "chr1\t$start\t$end\tC$v3\t$v1\t+\n";
      }
    }
    else
    {
      print $outfile "chr1\t$start\t$end\tC$v3\t$v1\t+\n";
    }
  }
  close PARSERNASEQ;
}
elsif ($cmd eq "combine")
{
  # Find the number of parsernaseq files.
  my $yes = 1;
  my $i = 0;
  while ($yes == 1)
  {
    $i++;
    my $f = $parsernaseq.$i;
    unless (-e $f)
    {
      $yes = 0;
    }
  }
  print $i, "\n";
  my $numberParseRNAseqFile = $i - 1;

  for ($i = 1; $i <= $numberParseRNAseqFile; $i++)
  {
    my $f = $parsernaseq.$i;

    my $topLevel = 1;
    open PARSERNASEQ, $f or die "cannot open < $f $!";
    while (<PARSERNASEQ>)
    {
      chomp;
      my @e = split /\t/;
      my $v2    = $e[8];
      $v2 =~ /Name=C(\d+)\_/; 
      my $v3 = $1;
      if ($topLevel < $v3)
      {
        $topLevel = $v3;
      }
    }
    close PARSERNASEQ;
    
    open PARSERNASEQ, $f or die "cannot open < $f $!";
    while (<PARSERNASEQ>)
    {
      chomp;
      my @e = split /\t/;
      my $start = $e[3];
      my $end   = $e[4];
      my $v1    = $e[5];
      my $v2    = $e[8];
      $v2 =~ /Name=C(\d+)\_/; 
      my $v3 = $1;
      if ($i < $numberParseRNAseqFile)
      {
        if ($v3 == $topLevel)
        {
          print $outfile "chr1\t$start\t$end\tC$i\t$v1\t+\n";
        }
      }
      else
      {
        if ($v3 > 4)
        {
          print $outfile "chr1\t$start\t$end\tC$i-$v3\t$v1\t+\n";
        }
      }
    }
    close PARSERNASEQ;
  }
}
elsif ($cmd eq "remove")
{
  # Find the highest level.
  my $topLevel = 1;
  open PARSERNASEQ, $parsernaseq or die "cannot open < $parsernaseq $!";
  while (<PARSERNASEQ>)
  {
    chomp;
    my @e = split /\t/;
    my $v2    = $e[8];
    $v2 =~ /Name=C(\d+)\_/; 
    my $v3 = $1;
    if ($topLevel < $v3)
    {
      $topLevel = $v3;
    }
  }
  close PARSERNASEQ;

  # Find the regions with the highest level, and 
  # zero the pileup values in the regions
  my $iPileup = 1;
  open PILEUP, $params{inpileup} or die "cannot open < $params{pileup} $!";
  open PARSERNASEQ, $parsernaseq or die "cannot open < $parsernaseq $!";
  while (<PARSERNASEQ>)
  {
    chomp;
    my @e = split /\t/;
    my $start = $e[3];
    my $end   = $e[4];
    my $v1    = $e[5];
    my $v2    = $e[8];
    $v2 =~ /Name=C(\d+)\_/; 
    my $v3 = $1;
    if ($v3 == $topLevel)
    {
      my $l;
      for (; $iPileup < $start; $iPileup++)
      {
        $l = <PILEUP>;
        print $outfile $l;
      }
      for (; $iPileup <= $end; $iPileup++)
      {
        $l = <PILEUP>;
        $l =~ /^(\d+)\s+/;
        print $outfile "$1\t0\n"; 
      }
    }
  }
  while (my $l = <PILEUP>) {
    print $outfile $l;
  }
  close PARSERNASEQ;
  close PILEUP;

  # Use the new pileup file to remove genes in the feature file.
  my @pileupValue;
  push @pileupValue, 0; # The 1-index pileupValue
  if (exists $params{out})
  {
    close $outfile;
  }
  $out = $params{out};
  open NEWPILEUP, $out or die "cannot open < $out: $!";
  while (<NEWPILEUP>)
  {
    /^(\d+)\s+(\d+)/;
    push @pileupValue, $2;
  }
  close NEWPILEUP;
  # SMU.01    0    193    1552    1360
  $feature =~ /(.+feature-genome\.out-parsernaseq)(\d+)/;
  my $iterationID = $2 + 1;
  my $newFeature = "$1$iterationID";
  open NEWFEATURE, ">", $newFeature or die "cannot open > $newFeature: $!";
  open FEATURE, $feature or die "cannot open < $feature: $!";
  while (my $l = <FEATURE>)
  {
    chomp $l;
    my @e = split /\s+/, $l;
    my $start = $e[2] + 1;
    my $end = $e[3];
    my $length = $end - $start + 1;
    my $numberZero = 0;
    for (my $i = $start; $i <= $end; $i++)
    {
      if ($pileupValue[$i] == 0)
      {
        $numberZero++; 
      }
    }
    if ($numberZero/$length < 0.1)
    {
      print NEWFEATURE "$l\n";
    }
    else
    {
      print STDERR "$l\t$numberZero/$length\n";
    }
  }
  close FEATURE;
  close NEWFEATURE;
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
    # print "$t->{name}: ";
    my @genesT;
    if (exists $t->{gene})
    {
      @genesT = @{$t->{gene}};
      for (my $j = 0; $j <= $#genesT; $j++)
      {
        my $g = $genesT[$j];
        # print ":$g->{name}";
      } 
    }
    # print "\n";
  }
  # print "\n++++++++++++++++++++++++++++++++\n";

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
    # print "\n";
  }

  for (my $i = 0; $i <= $#m2; $i++)
  {
    my $t = $m2[$i];
    print $outfile "chr1\t$t->{start}\t$t->{end}\t$t->{name}\t$t->{score}\t$t->{strand}\n";
    # print $outfile "chr1\t$t->{start}\t$t->{end}\t$t->{name}\t$t->{score}\t$t->{strand}\n";
  }
}
elsif ($cmd eq "adjust")
{
  my @e = rnaseqWiggleParse ($end);
  my @o = rnaseqBedParse ($operon);

  my $genomeLength = $#e + 1;

  # Use $peakcutoff and $sizecutoff
  # Use $windowsize
  for (my $i = 0; $i <= $#o; $i++)
  {
    my $t = $o[$i];
    my @peaks;
    if ($t->{strand} eq '+')
    {
      my $walkstart = $t->{start} - $windowsize;
      my $walkend = $t->{start} + $windowsize;
      if ($walkstart < 0)
      {
        $walkstart = 0;
      }
      if ($walkend >= $genomeLength)
      {
        $walkend = $#e;
      }
      for (my $j = $walkstart; $j < $walkend; $j++)
      {
        if ($e[$j] > $peakcutoff)
        {
          # Check if there are no sign change for the next sizecutoff positions.
          my $signchange = 0;
          my $jend = $j + $sizecutoff;
          $jend = $#e if $jend >= $genomeLength;
          for (my $k = $j; $k < $jend; $k++)
          {
            if ($e[$k] < 0)
            {
              $signchange = 1;
            }
          }
          if ($signchange == 0)
          {
            push @peaks, $j;
          }
        }
      }
      my @peakFound = @peaks;
      $t->{peak} = \@peakFound;
    }
    else
    {
      # Note that this is the same in the code above.
      # Do this for the negative strand.
      my $walkstart = $t->{end} - $windowsize;
      my $walkend = $t->{end} + $windowsize;
      if ($walkstart < 0)
      {
        $walkstart = 0;
      }
      if ($walkend >= $genomeLength)
      {
        $walkend = $#e;
      }
      for (my $j = $walkend; $j > $walkstart; $j--)
      {
        if ($e[$j] < 0 - $peakcutoff)
        {
          # Check if there are no sign change for the next sizecutoff positions.
          my $signchange = 0;
          my $jstart = $j - $sizecutoff;
          $jstart = 0 if $jstart < 0;
          for (my $k = $j; $k > $jstart; $k--)
          {
            if ($e[$k] > 0)
            {
              $signchange = 1;
            }
          }
          if ($signchange == 0)
          {
            push @peaks, $j;
          }
        }
      }
      my @peakFound = @peaks;
      $t->{peak} = \@peakFound;
    }
    # print STDERR "Processed ($i/$#o)\r";
  }

  for (my $i = 0; $i <= $#o; $i++)
  {
    my $t = $o[$i];
    my @peaks = @{$t->{peak}};
    if ($#peaks >= 0)
    {
      if ($t->{strand} eq '+')
      {
        $t->{thickstart} = $peaks[0];
        $t->{thickend} = $t->{end};
        if ($t->{thickstart} > $t->{thickend})
        {
          $t->{thickstart} = $t->{start};
        }
      }
      else
      {
        $t->{thickstart} = $t->{start};
        $t->{thickend} = $peaks[$#peaks];
        if ($t->{thickstart} > $t->{thickend})
        {
          $t->{thickend} = $t->{end};
        }
      }
    }
    else
    {
      $t->{thickstart} = $t->{start};
      $t->{thickend} = $t->{end};
    }

    print STDERR "$t->{name}\t$#peaks\n";
    #print $outfile "chr1\t$t->{start}\t$t->{end}\t$t->{name}\t$t->{score}\t$t->{strand}\t$t->{thickstart}\t$t->{thickend}\t255,0,0\n";
    print $outfile "chr1\t$t->{thickstart}\t$t->{thickend}\t$t->{name}\t$t->{score}\t$t->{strand}\n";
    #for (my $j = 0; $j <= $#peaks; $j++)
    #{
      #print $outfile "$peaks[$j]\t";
    #}
    #print $outfile "\n";
  }
}
elsif ($cmd eq "slope")
{
  my $windowsize2 = 10;
  my $nonzerosize = 1;
  my @e = rnaseqWiggleParse ($end);
  my @o = rnaseqBedParse ($operon);

  my $genomeLength = $#e + 1;
  my @map1 = (0) x $genomeLength;
  my @win1 = (0) x $windowsize;

  ################################################
  # The first windowsize-many elements are stored.
  for (my $i = 0; $i < $windowsize; $i++)
  {
    $win1[$i] = $e[$i]; 
  }

  ################################################
  # Find the sign of the windowsize-many elements.
  my %numSign;
  my $sign;
  $numSign{plus} = 0;
  $numSign{mnus} = 0;
  $numSign{zero} = 0;
  for (my $i = 1; $i < $windowsize; $i++)
  {
    if ($win1[$i] > 0)
    {
      $numSign{plus}++;
    }
    elsif ($win1[$i] < 0)
    {
      $numSign{mnus}++;
    }
    else 
    {
      $numSign{zero}++;
    }
  }

  ################################################
  # Start to check the windowsize-many elements.
  for (my $i = 0; $i < $genomeLength - $windowsize - 1; $i++)
  {
    if ($numSign{mnus} == 0)
    {
      my $endrange = $windowsize2 - 1;
      my @win2 = @win1[0..$endrange];
      my $numberNonzero = arrayNumberNonzero (@win2);
      my $slopeValue = 0;
      if ($numberNonzero > $nonzerosize)
      {
        $slopeValue = arrayNumberSlope (@win2);
      }
      $map1[$i] = $slopeValue;
    }
    elsif ($numSign{plus} == 0)
    {
      my @win1r = arrayNumberNegate(reverse (@win1));
      my $endrange = $windowsize2 - 1;
      my @win2 = @win1r[0..$endrange];
      my $numberNonzero = arrayNumberNonzero (@win2);
      my $slopeValue = 0;
      if ($numberNonzero > $nonzerosize)
      {
        $slopeValue = arrayNumberSlope (@win2);
      }
      $map1[$i+$windowsize-1] = -$slopeValue;
    }

    ###############################################
    # Pop and push elements.
    # Remove the first element.
    # print STDERR "$i\n";
    # last if ($i > 100);
    if ($win1[0] > 0)
    {
      $numSign{plus}--;
    }
    elsif ($win1[0] < 0)
    {
      $numSign{mnus}--;
    }
    else 
    {
      $numSign{zero}--;
    }
    shift @win1;
    push @win1, $e[$i+$windowsize];
    if ($win1[$#win1] > 0)
    {
      $numSign{plus}++;
    }
    elsif ($win1[$#win1] < 0)
    {
      $numSign{mnus}++;
    }
    else 
    {
      $numSign{zero}++;
    }
  }

  ################################################
  # Find the position at which the slope is the largest
  # among the sites within 200 bp up- and down-stream 
  # of the start and end of transcripts.
  for (my $i = 0; $i <= $#o; $i++)
  {
    my $t = $o[$i];
    my @peaks;
    my @valleys;

    # Start of an operon
    my $walkstart = ($t->{start} - 1) - $updowncutoff;
    my $walkend = ($t->{start} - 1) + $updowncutoff;
    $walkstart = 0 unless $walkstart >= 0;
    $walkend = $genomeLength - 1 unless $walkend < $genomeLength;
    my @win3 = @map1[$walkstart..$walkend];
    my $idxMax = 0;
    $win3[$idxMax] > $win3[$_] or $idxMax = $_ for 1 .. $#win3;
    @peaks = ();
    if ($win3[$idxMax] > 0)
    {
      push @peaks, $walkstart + $idxMax; 
    }
    my @peakFound = @peaks;
    $t->{peak} = \@peakFound;

    # End of an operon
    $walkstart = ($t->{end} - 1) - $updowncutoff;
    $walkend = ($t->{end} - 1) + $updowncutoff;
    $walkstart = 0 unless $walkstart >= 0;
    $walkend = $genomeLength - 1 unless $walkend < $genomeLength;
    @win3 = @map1[$walkstart..$walkend];
    $idxMax = 0;
    $win3[$idxMax] < $win3[$_] or $idxMax = $_ for 1 .. $#win3;
    @valleys = ();
    if ($win3[$idxMax] < 0)
    {
      push @valleys, $walkstart + $idxMax; 
    }
    my @valleyFound = @valleys;
    $t->{valley} = \@valleyFound;

  }
  
  for (my $i = 0; $i <= $#map1; $i++)
  {
    print $map1[$i], "\n";
  } 

  for (my $i = 0; $i <= $#o; $i++)
  {
    my $t = $o[$i];
    my @peaks = @{$t->{peak}};
    my @valleys = @{$t->{valley}};
    if ($#peaks >= 0)
    {
      $t->{thickstart} = $peaks[0];
    }
    else
    {
      $t->{thickstart} = $t->{start};
    }
    if ($#valleys >= 0)
    {
      $t->{thickend} = $valleys[0] + 1;
    }
    else
    {
      $t->{thickend} = $t->{end};
    }
    unless ($t->{thickstart} < $t->{thickend})
    {
      $t->{thickstart} = $t->{start};
      $t->{thickend} = $t->{end};
    }

    print STDERR "$t->{name}\t$#peaks\t$#valleys\n";
    #print $outfile "chr1\t$t->{start}\t$t->{end}\t$t->{name}\t$t->{score}\t$t->{strand}\t$t->{thickstart}\t$t->{thickend}\t255,0,0\n";
    print $outfile "chr1\t$t->{thickstart}\t$t->{thickend}\t$t->{name}\t$t->{score}\t$t->{strand}\n";
    #for (my $j = 0; $j <= $#peaks; $j++)
    #{
      #print $outfile "$peaks[$j]\t";
    #}
    #print $outfile "\n";
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

perl pl/transcript-parsernaseq.pl combine -parsernaseq out.parsernaseq -out bed.file

perl pl/transcript-parsernaseq.pl operon -feature feature-genome.out-geneonly -parsernaseq out.parsernaseq -out operon.file

perl pl/transcript-parsernaseq.pl adjust -end FASTQ01-end.wig -operon operon.file -out adjusted.file

perl pl/transcript-parsernaseq.pl slope -end FASTQ01-end.wig -operon operon.file -out adjusted.file

perl pl/transcript-genecoverage.pl remove -feature feature.file -inpileup FASTQNUM.pileup -parsernaseq FASTQNUM.parsernaseq 

=head1 DESCRIPTION

transcript-parsernaseq will help you to pre- and post-process files of
ParseRNASeq.

Command:
  pileup - Using FASTQ01.wig file that is created by samtools mpileup to create
  a pileup file for ParaseRNAseq.

  gene - Using feature-genome.out-geneonly file that is created by menu
  feature-genome I create a gene annotation file for ParseRNASeq.

  bed - Using ParseRNAseq output to create a bed file.

  gff - Using ParseRNAseq gff output to create a bed file.

  combine - Using ParseRNAseq gff output files to create one bed file.

  operon - Using ParseRNAseq gff output to create an operon file.

  adjust - Operons are merely groups of genes. Their start positions are those
  of the first genes in the operons. RNA-Seq data ends information can be useful
  in adjusting the start positions. Menu bwa-pos2wig should be used to have a
  FASTQ01-end.wig file.

  slope - For each position of a genome I compute the slope of increasing or
  decreasing peaks using FASTQ01-end.wig file. A site with this pattern is a
  candidate transcription start or end site. A site should pass filters: 1) the
  following 95 base pairs must be the same sign including the site itself. 2)
  Five of the the following 10 bases must be nonzero. 3) Using the following 10
  bases I compute the slope.

  remove - Remove regions of pileups where the highest level of expressed
  transcripts are predicted. Or, make all values zeros at the regions of predicted
  transcripts with highest expression level category. ParseRNAseq -force_gp
  option uses annotated genes as transcribed regions. I do not want genes with
  no pileup values to be used as transcribed ones. This command creates two
  files: one is a new feature file, and another is a new pileup file.

=head1 OPTIONS

=over 8

=item B<-feature> <file>

Features of a genome. 

=item B<-wiggle> <file>

The wiggle file of coverage values.

=item B<-parsernaseq> <file>

The output file from ParseRNAseq.

=item B<-end> <file>

The output file from bwa-pos2wig's end command.

=item B<-operon> <file>

The output file from operon command.

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
