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

###############################################################################
# COMMAND LINE
###############################################################################
$| = 1; # Do not buffer output
my $VERSION = 'de-count.pl 1.0';

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
            'first',
            'singlegenome',
            'shortread=s',
            'genepos=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $shortread;
my $shortreadfile;
my $genepos;
my $geneposfile;
my $out;
my $outfile;

if (exists $params{out})
{
  $out = $params{out};
  open ($outfile, ">", $out) or die "cannot open > $out: $!";
}
else
{
  $outfile = *STDOUT;   
}

if (exists $params{shortread})
{
  $shortread = $params{shortread};
}
else
{
  &printError("short read file is missing");
}

if (exists $params{genepos})
{
  $genepos = $params{genepos};
}
else
{
  &printError("genepos file is missing");
}

###############################################################################
# DATA PROCESSING
###############################################################################

if ($cmd eq "overlap")
{
  # my @shortreads = rnaseqPosParse ($shortread); This could be slow if the
  # number of short reads exceeds millions.
  my @genes = rnaseqPosParse ($genepos);
  for (my $i = 0; $i < $#genes; $i++)
  {
    my $gi = $genes[$i];
    for (my $j = $i + 1; $j <= $#genes; $j++)
    {
      # genei :     |---------------|
      # case 1:                 <------> 
      # case 2:  <------> 
      # case 3:        <------>
      # case 4:  <--------------------->
      my $gj = $genes[$j];
      if ($gi->{chr} eq $gj->{chr})
      {
        my $v = 0;
        if ($gi->{start} <= $gj->{start} and $gj->{start} <= $gi->{end})
        { 
          $v++; # case 1 or 3
        }
        elsif ($gi->{start} <= $gj->{end} and $gj->{end} <= $gi->{end})
        { 
          $v++; # case 2 or 3
        }
        elsif ($gj->{start} <= $gi->{start} and $gi->{end} <= $gj->{end})
        { 
          $v++; # case 4
        }
        else
        {
          # die "Impossible case: s-start $s->{start}, s-end $s->{end}, g-start $g->{start}, g-end  $g->{end}"
        }
        if ($v > 0)
        {
          my $distance;
          if ($gi->{start} < $gj->{start})
            {
              $distance = $gi->{end} - $gj->{start};
            }
          else
            {
              $distance = $gj->{end} - $gi->{start};
            }
          print STDERR "$gi->{name} and $gj->{name} are overlapped:";
          print STDERR "$gi->{chr} $gi->{start} $gi->{end}: "; 
          print STDERR "$gj->{chr} $gj->{start} $gj->{end}: $distance\n"; 
        }
      }
    }
  }
}
elsif ($cmd eq "join")
{
  # my @shortreads = rnaseqPosParse ($shortread); This could be slow if the
  # number of short reads exceeds millions.
  my @genes = rnaseqPosParse ($genepos);

  my $c = 0;
  open SHORTREAD, $shortread or die "cannot open < $shortread";
  my $line;
  if (exists $params{first})
  {
    $line = <SHORTREAD>;
  }
  while (<SHORTREAD>)
  {
    chomp;
    my @e = split /\t/;
    my $s = {};
    $s->{name}  = $e[0]; 
    $s->{chr}   = $e[1]; 
    $s->{start} = $e[2]; 
    $s->{end}   = $e[3]; 
    # gene  :     |---------------|
    # case 1:                 <------> 
    # case 2:  <------> 
    # case 3:        <------>
    # case 4:  <--------------------->
    unless ($s->{start} < $s->{end})
    {
      die "short read start must be less than its end";
    }

    for (my $i = 0; $i <= $#genes; $i++)
    {
      my $g = $genes[$i];

      # When UA159 is used as a reference genome, there is no need to check if
      # the chromosome of a short read and that of a gene are the same.
      # When OMZ175 or SMU86 genome is used as a reference genome, we have to do
      # this because the genome consists of contigs. Effectively, there are
      # multiple chromosomes.
      # We have to change feature-genome.out-geneonly or
      # feature-genome.pl script.
      if (exists $params{singlegenome} || $g->{chr} eq $s->{chr}) # Why did I comment this?
      {
        my $v = $g->{count};
        if ($g->{start} <= $s->{start} and $s->{start} <= $g->{end})
        { 
          $v++; # case 1 or 3
        }
        elsif ($g->{start} <= $s->{end} and $s->{end} <= $g->{end})
        { 
          $v++; # case 2 or 3
        }
        elsif ($s->{start} <= $g->{start} and $g->{end} <= $s->{end})
        { 
          $v++; # case 4
        }
        else
        {
          # die "Impossible case: s-start $s->{start}, s-end $s->{end}, g-start $g->{start}, g-end  $g->{end}"
        }
        $g->{count} = $v;
      }
    }
    $c++;
    if ($c % 1000000 == 0)
    #if ($c % 100000 == 0)
    {
      print STDERR "Read $c\r";
    }
  }
  for (my $i = 0; $i <= $#genes; $i++)
  {
    my $g = $genes[$i];
    # print $outfile "$g->{name}\t$g->{count}\n";
    print $outfile "$g->{count}\n";
  }
}
elsif ($cmd eq "strict")
{
  # my @shortreads = rnaseqPosParse ($shortread); This could be slow if the
  # number of short reads exceeds millions.
  my @genes = rnaseqPosParse ($genepos);

  my $c = 0;
  open SHORTREAD, $shortread or die "cannot open < $shortread";
  while (<SHORTREAD>)
  {
    chomp;
    my @e = split /\t/;
    my $s = {};
    $s->{name}  = $e[0]; 
    $s->{chr}   = $e[1]; 
    $s->{start} = $e[2]; 
    $s->{end}   = $e[3]; 
    # gene  :     |---------------|
    # case 1:                 <------> 
    # case 2:  <------> 
    # case 3:        <------>
    # case 4:  <--------------------->
    unless ($s->{start} < $s->{end})
    {
      die "short read start must be less than its end";
    }

    for (my $i = 0; $i <= $#genes; $i++)
    {
      my $g = $genes[$i];
      my $v = $g->{count};
      if ($g->{chr} eq $s->{chr})
      {
        if ($g->{start} <= $s->{start} and $s->{start} <= $g->{end}
            and $g->{start} <= $s->{end} and $s->{end} <= $g->{end})
        { 
          $v++; # only case 3
        }
        else
        {
          # die "Impossible case: s-start $s->{start}, s-end $s->{end}, g-start $g->{start}, g-end  $g->{end}"
        }
      }
      $g->{count} = $v;
    }
    $c++;
    if ($c % 1000000 == 0)
    {
      print STDERR "$c\r";
    }
  }
  for (my $i = 0; $i <= $#genes; $i++)
  {
    my $g = $genes[$i];
    print $outfile "$g->{name}\t$g->{count}\n";
  }
}

if (exists $params{out})
{
  close $outfile;
}
__END__
=head1 NAME

de-count - count of short reads on genes

=head1 VERSION

de-count 1.0

=head1 SYNOPSIS

perl de-count.pl join -shortread FASTQ01-sum.pos -genepos gene.pos

=head1 DESCRIPTION

de-count will help you to count short reads mapped on genes.

Command:
  join    - count short reads mapped on genes using shortreadfile and geneposfile.
  strict  - count short reads mapped only within genes using shortreadfile and geneposfile.
  overlap - shows what genes are overlapped

  shortreadfile should contain 4 columns: short read name, chromosome, start,
  and end positions.

  geneposfile should contain 4 columns: gene name, chromosome, start, and end
  positions.

=head1 OPTIONS

=over 8

=item B<-shortread> <file>

A file should should contain 4 columns: short read name, chromosome, start, and
end positions.

=item B<-genepos> <file>

A file should contain 4 columns: gene name, chromosome, start, and end
positions.

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
