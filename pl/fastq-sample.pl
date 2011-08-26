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
my $VERSION = 'fastq-sample.pl 1.0';

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
            'fastq=s',
            'out=s',
            'interval=i',
            'proportion=i',
            'samplesize=i',
            'cutsize=i',
            'keepsize=i',
            'removesize=i',
            'maxsize=i',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $fastq;
my $interval;
my $proportion;
my $samplesize;
my $out;
my $outfile;
my $cutsize;
my $keepsize;
my $removesize;
my $maxsize;

if (exists $params{fastq}) 
{
  $fastq = $params{fastq};
}
else
{
  &printError("fastq file is missing");
}

if (exists $params{samplesize}) 
{
  $samplesize = $params{samplesize};
  $proportion = $samplesize / 100;
}

if (exists $params{proportion}) 
{
  $proportion = $params{proportion};
  $proportion /= 100;
}

if (exists $params{interval}) 
{
  $interval = $params{interval};
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

if (exists $params{keepsize}) 
{
  $keepsize = $params{keepsize};
}

if (exists $params{maxsize}) 
{
  $maxsize = $params{maxsize};
}

if (exists $params{removesize}) 
{
  $removesize = $params{removesize};
}

if (exists $params{cutsize}) 
{
  $cutsize = $params{cutsize};
}

if ($cmd eq "sample")
{
  unless (exists $params{interval} or exists $params{proportion}
          or exists $params{samplesize}) 
  {
    &printError("command sample needs -samplesize or -interval or -proportion");
  }
}
elsif ($cmd eq "cut")
{
  unless (exists $params{cutsize}) 
  {
    &printError("command cut needs -cutsize");
  }
}
elsif ($cmd eq "keep")
{
  unless (exists $params{keepsize}) 
  {
    &printError("command keep needs -keepsize");
  }
}
elsif ($cmd eq "remove")
{
  unless (exists $params{removesize}) 
  {
    &printError("command remove needs -removesize");
  }
}
elsif ($cmd eq "max")
{
  unless (exists $params{maxsize}) 
  {
    &printError("command max needs -maxsize");
  }
}

###############################################################################
# DATA PROCESSING
###############################################################################
if ($cmd eq "sample")
{
  my $skip; 
  # $random = int( rand( $Y-$X+1 ) ) + $X;
  my $X;
  my $Y;
  my $Z;
  if (defined $interval)
  {
    $X = 1;
    $Y = $interval*2;
  }
  elsif (defined $proportion) # See pl/proportion.pl
  {
    $X = int(1/$proportion);
    $Y = 1/$proportion - $X;
    $Z = $Y/$X;
  }

  my $line;
  my $line1;
  my $line2;
  my $line3;
  my $line4;
  open FASTQ, "zcat $fastq|" or die "cannot open < $fastq $!";
  while ($line1 = <FASTQ>)
  {
    $line2 = <FASTQ>;
    $line3 = <FASTQ>;
    $line4 = <FASTQ>;
    if (defined $interval)
    {
      my $skip = int (rand($Y-$X+1)) + $X;
    }
    elsif (defined $proportion)
    {
      $skip = 2*$X;
      my $r1 = rand();
      if ($r1 >= $Z)
      {
        my $r2 = int (rand(2*$X-1)) + 1;
        $skip = $r2;
      }
    }
    for (my $i = 0; $i < $skip - 1; $i++)
    # for (my $i = 1; $i < $skip; $i++)
    {
      $line1 = <FASTQ>;
      unless ($line1)
      {
        last;
      }
      $line2 = <FASTQ>;
      $line3 = <FASTQ>;
      $line4 = <FASTQ>;
    }
    if ($line1)
    {
      print $outfile $line1;
      print $outfile $line2;
      print $outfile $line3;
      print $outfile $line4;
    }
    
    #for my $j (1..4) 
    #{
    #  $line = <FASTQ>; 
    #  unless ($line)
    #  {
    #    last;
    #  }
    #  print $outfile $line;
    #}
  }
  close FASTQ;
}
elsif ($cmd eq "cut")
{
  my $s = 0;
  my ($line1, $line2, $line3, $line4);
  open FASTQ, "zcat $fastq|" or die "cannot open < $fastq $!";
  while ($line1 = <FASTQ>)
  {
    $s++;
    $line2 = <FASTQ>;
    chomp $line2;
    my $l = length $line2;
    my $shortenedLine2; 
    $line3 = <FASTQ>;
    $line4 = <FASTQ>;
    chomp $line4;
    my $shortenedLine4;
    $shortenedLine2 = substr ($line2, 0, $l - $cutsize);
    $shortenedLine4 = substr ($line4, 0, $l - $cutsize);
    print $outfile $line1;
    print $outfile "$shortenedLine2\n";
    print $outfile $line3;
    print $outfile "$shortenedLine4\n";
    if ($s % 1000000 == 0)
    {
      print STDERR "Read $s\r";
    }
  }
  print STDERR "Finished Read $s for cutting $cutsize tailing bases.\n";
}
elsif ($cmd eq "max")
{
  my $s = 0;
  my ($line1, $line2, $line3, $line4);
  open FASTQ, "zcat $fastq|" or die "cannot open < $fastq $!";
  while ($line1 = <FASTQ>)
  {
    $s++;
    $line2 = <FASTQ>;
    chomp $line2;
    my $l = length $line2;
    my $shortenedLine2; 
    $line3 = <FASTQ>;
    $line4 = <FASTQ>;
    chomp $line4;
    my $shortenedLine4;
    if ($l > $maxsize)
    {
      $shortenedLine2 = substr ($line2, 0, $maxsize);
      $shortenedLine4 = substr ($line4, 0, $maxsize);
    }
    else
    {
      $shortenedLine2 = $line2;
      $shortenedLine4 = $line4;
    }
    print $outfile $line1;
    print $outfile "$shortenedLine2\n";
    print $outfile $line3;
    print $outfile "$shortenedLine4\n";
    if ($s % 1000000 == 0)
    {
      print STDERR "Read $s\r";
    }
  }
  print STDERR "Finished Read $s for reads upto $maxsize size\n";
}
elsif ($cmd eq "keep")
{
  my $s = 0;
  my ($line1, $line2, $line3, $line4);
  open FASTQ, "zcat $fastq|" or die "cannot open < $fastq $!";
  while ($line1 = <FASTQ>)
  {
    $s++;
    $line2 = <FASTQ>;
    chomp $line2;
    my $l = length $line2;
    $line3 = <FASTQ>;
    $line4 = <FASTQ>;
    if ($l >= $keepsize)
    {
      print $outfile $line1;
      print $outfile "$line2\n";
      print $outfile $line3;
      print $outfile $line4;
    }
    if ($s % 1000000 == 0)
    {
      print STDERR "Read $s\r";
    }
  }
  print STDERR "Finished Read $s for keeping $keepsize bases.\n";
}
elsif ($cmd eq "remove")
{
  my $s = 0;
  my ($line1, $line2, $line3, $line4);
  open FASTQ, "zcat $fastq|" or die "cannot open < $fastq $!";
  while ($line1 = <FASTQ>)
  {
    $s++;
    $line2 = <FASTQ>;
    chomp $line2;
    my $l = length $line2;
    $line3 = <FASTQ>;
    $line4 = <FASTQ>;
    if ($l < $removesize)
    {
      print $outfile $line1;
      print $outfile "$line2\n";
      print $outfile $line3;
      print $outfile $line4;
    }
    if ($s % 1000000 == 0)
    {
      print STDERR "Read $s\r";
    }
  }
  print STDERR "Finished Read $s for removing minimum length of $removesize bases.\n";
}

if (exists $params{out})
{
  close $outfile;
}
__END__
=head1 NAME

fastq-sample - Subsample a FASTQ file

=head1 VERSION

fastq-sample 1.0

=head1 SYNOPSIS

perl fastq-sample.pl sample -fastq file -out file -interval number
perl fastq-sample.pl sample -fastq file -out file -proportion number
perl fastq-sample.pl cut -fastq file -out file -cutsize 10
perl fastq-sample.pl keep -fastq file -out file -keepsize 10
perl fastq-sample.pl remove -fastq file -out file -removesize 10

=head1 DESCRIPTION

fastq-sample will help you to subsample a FASTQ file to create a smaller FASTQ
file from the original one.

=head1 OPTIONS

command: sample, cut
  sample - sample short reads
  cut    - truncate tailing nucleotides from short reads
  keep   - keep only short reads with minimum length
  remove - remove reads greater than or equal to a given length
  max    - keep short reads with maximum length

=over 8

=item B<-fastq> <file>

A FASTQ file name.

=item B<-removesize> <number>

Number of minimum lengths of reads that are removed.

=item B<-maxsize> <number>

Number of maximum lengths of reads to keep.

=item B<-cutsize> <number>

Number of tailing nucleotides.

=item B<-keepsize> <number>

Minimum number of nucleotides of short reads

=item B<-proportion> <number>

A number between 0 and 100.

=item B<-interval> <number>

Average interval of short reads to jump for next samples.

=item B<-out> <file>

An output file name.

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


