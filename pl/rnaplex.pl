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
my $VERSION = 'rnaplex.pl 1.0';

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
            'sortbyenergy',
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

################################################################################
## DATA PROCESSING
################################################################################
if ($cmd eq "rank")
{
  # >target
  # >query
  # (((.(((((((......((((((((&)))))))).....))))))).))) 106,130 :  13,36  (-4.49 = -21.00 +  7.46 +  9.05)
  my $line;
  my %queryTarget;
  while ($line = <$infile>)
  {
    chomp $line;
    $line = substr $line, 1;
    my $target = $line;
    $line = <$infile>;
    chomp $line;
    $line = substr $line, 1;
    my $query = $line;
    $line = <$infile>;
    chomp $line;
    $line =~ /\((-\d+\.\d+)\s+=/;
    my $value = $1;

    if (exists $queryTarget{$query})
    {
      my $rec = $queryTarget{$query};
      if (exists $rec->{$target})
      {
        die "One target for a query must exist.";
      }
      $rec->{$target} = $value;
    }
    else
    {
      my $rec = {};
      $queryTarget{$query} = $rec;
      $rec->{$target} = $value;
    }
  }

  # Sort by values.
  if (exists $params{sortbyenergy})
  {
    foreach my $query (keys %queryTarget)
    {
      my %targets = %{$queryTarget{$query}};
      my $i = 0;
      print $outfile "$query";
      foreach my $value (sort {$targets{$a} <=> $targets{$b} } keys %targets)
      {
        $i++;
        print $outfile "\t$value:$targets{$value}";
        if ($i > 5)
        {
          # last;
        }
      }
      print $outfile "\n";
    }
  }
  else
  {
    my @queries = sort {$queryTarget{$b} cmp $queryTarget{$a}} keys (%queryTarget);
    my $firstQuery = $queries[0];
    my %targets = %{$queryTarget{$firstQuery}};

    print $outfile "target";
    foreach my $q (@queries)
    {
      print $outfile "\t$q";
    }
    print $outfile "\n";

    foreach my $t (keys %targets)
    {
      print $outfile "$t";
      foreach my $q (@queries)
      {
        print $outfile "\t$queryTarget{$q}{$t}";
      }
      print $outfile "\n";
    }
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

rnaplex.pl - Parse RNAplex output files.

=head1 VERSION

rnaplex.pl 1.0

=head1 SYNOPSIS

perl rnaplex.pl rank -in target.rnaplex -out target.rnaplex.sorted

=head1 DESCRIPTION

rnaplex.pl parses output files of RNAplex.

=head1 OPTIONS

  command: 

  rank - Read an output file of RNAplex to list targets in order for each query
  small RNA.

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

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

=head1 BUGS

If you find a bug please post a message rnaseq_analysis project at codaset dot
com repository so that I can make rnaplex.pl better.

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
