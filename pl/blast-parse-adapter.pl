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

###############################################################################
# COMMAND LINE
###############################################################################
$| = 1; # Do not buffer output
my $VERSION = 'blast-parse-adapter.pl 1.0';

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
            'in=s',
            'out=s',
            'reads=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $out;
my $in;
my $outfile;
my $infile;
my $reads;

if (exists $params{reads})
{
  $reads = $params{reads};
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

if (exists $params{in})
{
  $in = $params{in};
  open ($infile, "<", $in) or die "cannot open > $in: $!";
}
else
{
  $infile = *STDIN;   
}

if ($cmd eq "filterfastq")
{
  unless (exists $params{reads})
  {
    &printError("reads is missing");
  }
}

###############################################################################
# DATA PROCESSING
###############################################################################

if ($cmd eq "plusplus")
{
  my $line;
  # Go to the first alignment.
  while ($line = <$infile>)
  {
    chomp $line;
    last if $line =~ /^>/;
  }

  # Start to find reads and adapter's position
  my %aln;
  $aln{readsName} = $line;
  while ($line = <$infile>)
  {
    chomp $line;
    if ($line =~ /^>/)
    {
      # Print the previous reads
      if ($aln{identities} == 100 and $aln{adapterPos} > 10)
      {
        print $outfile "$aln{readsName}\t$aln{adapterPos}\n";
      }

      # Start reading a new short read.
      $aln{readsName} = $line;
    }
    if ($line =~ /Identities\s+=\s+\d+\/\d+\s+\((\d+)\%\)/)
    {
      $aln{identities} = $1;
    }
    if ($line =~ /Strand\s+=\s+Plus\s+\/\s+Plus/)
    {
      $aln{strand} = "okay";
    }
    if ($line =~ /Query:\s+1\s+/)
    {
      $line = <$infile>;
      $line = <$infile>;
      if ($line =~ /Sbjct:\s+(\d+)\s+/)
      {
        $aln{adapterPos} = $1;
      }
      else
      {
        die "Error: There must be Sbjct";
      }
    }
  }
  # Print the previous reads
  if ($aln{identities} == 100 and $aln{adapterPos} > 10)
  {
    print $outfile "$aln{readsName}\t$aln{adapterPos}\n";
  }
}
elsif ($cmd eq "filterfastq")
{
  # Read short reads
  my %filterRead;
  open READS, $reads or die "cannot open $reads $!";
  while (<READS>)
  {
    chomp;
    my @e = split /\t/;
    $e[0] = substr $e[0], 1;
    $filterRead{$e[0]} = $e[1];
  }
  close READS;

  my $s = 0;
  my $line1;
  my $line2;
  my $line3;
  my $line4;
  while ($line1 = <$infile>)
  {
    $s++;
    $line2 = <$infile>;
    $line3 = <$infile>;
    $line4 = <$infile>;
    
    foreach my $r (keys %filterRead)
    {
      # print STDERR "Line1: $line1";
      # print STDERR "Line1-R: $r\n";
      my $ii = index ($line1, $r);
      if ($ii >= 0)
      {
        my $shortenLine2 = substr ($line2, 0, $filterRead{$r} - 1);
        my $shortenLine4 = substr ($line2, 0, $filterRead{$r} - 1);
        print $outfile $line1;
        print $outfile "$shortenLine2\n"; 
        print $outfile $line3;
        print $outfile "$shortenLine4\n"; 
        last;
      }
    }
    # print STDERR "\n\n";
    if ($s % 100 == 0)
    {
      print STDERR "Read $s\r";
    }
  }

  #foreach my $r (keys %filterRead)
  #{
    #print $outfile "$r\t$filterRead{$r}\n";
  #}
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

blast-parse-adapter - summary of BWA alignment

=head1 VERSION

blast-parse-adapter 1.0

=head1 SYNOPSIS

perl pl/blast-parse-adapter.pl plusplus -in blast.result -out reads.txt

zcat FASTQ.gz | perl pl/blast-parse-adapter.pl filterfastq -reads adapter.reads

=head1 DESCRIPTION

blast-parse-adapter will help you to parse BLAST results for positions at which
an adapter is aligned.

Command:
  plusplus - print short reads and adapter positions
  filterfastq - find short reads given by -reads option and shorten them

=head1 OPTIONS

=over 8

=item B<-in> <file>

If no input file is not specified, standard input is considered.

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
