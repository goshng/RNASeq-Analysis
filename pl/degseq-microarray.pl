#!/usr/bin/perl

#===============================================================================
#   Author: Sang Chul Choi, BSCB @ Cornell University, NY
#
#   File: degseq-microarray.pl
#   Date: Thu Jul  7 23:43:21 EDT 2011
#   Version: 1.0
#===============================================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

$| = 1; # Do not buffer output

my $VERSION = 'degseq-microarray.pl 1.0';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);        
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'degseq=s',
            'microarray=s',
            'reverse'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

degseq-microarray - join two files of degseq and microarray.

=head1 VERSION

degseq-microarray 1.0

=head1 SYNOPSIS

perl map2graph.pl [-degseq file] [-microarray file] 

=head1 DESCRIPTION

degseq-microarray will help you to merge degseq and microarray result.

=head1 OPTIONS

=over 8

=item B<-degseq> <file>

----
"GeneNames"	"value1"	"value2"	"log2(Fold_change)"	"log2(Fold_change) normalized"	"p-value"	"q-value(Benjamini et al. 1995)"	"q-value(Storey et al. 2003)"	"Signature(p-value < 0.001)"
SMU_06	6166	5635	0.129919279047591	0.343382760667889	0	0	0	TRUE
----
I use the values at the 4th and 5th columns. The gene name has a underscore at
the first column.

=item B<-microarray> <file>

----
0.0003325 43.729  SMU.1422  putative pyruvate dehydrogenase E1 component beta subunit)
----
I use the 2nd column and the 3rd one. Note that the gene name has a dot not
underscore.

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

require 'pl/sub-error.pl';

my $degseq;
my $microarray;
my $reverse = 0;

if (exists $params{degseq}) 
{
  $degseq = $params{degseq};
}
else
{
  &printError("degseq file is missing");
}

if (exists $params{microarray}) 
{
  $microarray = $params{microarray};
}
else
{
  &printError("microarray file is missing");
}

if (exists $params{reverse}) 
{
  $reverse = 1;
}

my %gene;
open DEGSEQ, $degseq or die "could not open $degseq $!";
while (<DEGSEQ>)
{
  my @e = split /\s+/;
  $gene{$e[0]}{DEGSEQ} = $e[4];
}
close DEGSEQ;

open MICRO, $microarray or die "could not open $microarray $!";
while (<MICRO>)
{
  my @e = split /\s+/;
  my $name = $e[2];
  $name =~ s/\./\_/;
  if ($reverse == 0)
  {
    $gene{$name}{MICRO} = log($e[1])/log(2);
  }
  else
  {
    $gene{$name}{MICRO} = -log($e[1])/log(2);
  }
}
close MICRO;

foreach my $f ( keys %gene ) {
  if (exists $gene{$f}{DEGSEQ} and exists $gene{$f}{MICRO})
  {
    print "$f\t$gene{$f}{DEGSEQ}\t$gene{$f}{MICRO}\n"
  }
}


