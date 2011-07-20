#!/usr/bin/perl
###############################################################################
# Copyright (C) 2011 Sang Chul Choi
#
# This file is part of RNASeq Analysis.
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

#===============================================================================
#   Author: Sang Chul Choi, BSCB @ Cornell University, NY
#
#   File: edgeR.pl
#   Date: Thu Jul  7 23:43:21 EDT 2011
#   Version: 1.0
#===============================================================================

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Scalar::Util qw(looks_like_number);
require 'pl/sub-error.pl';

$| = 1; # Do not buffer output

my $VERSION = 'edgeR.pl 1.0';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);        
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print $VERSION."\n"; exit; },
            'degseqout=s'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $degseqout;

if (exists $params{degseqout}) 
{
  $degseqout = $params{degseqout};
}
else
{
  &printError("degseqout directory is missing");
}


###############################################################################
# DATA PROCESSING
###############################################################################
sub raDegoutParse ($);
my $degoutput = raDegoutParse ("$degseqout/output_score.txt");

###############################################################################
# DATALST
###############################################################################
open DATALIST, ">", "$degseqout/datalist" or die "$!";
print DATALIST "files\tgroup\n";
print DATALIST "value1.txt\tvalue1\n";
print DATALIST "value2.txt\tvalue2\n";
close DATALIST;

for (my $i = 1; $i <= 2; $i++)
{
  open VALUE, ">", "$degseqout/value$i.txt" or die $!;
	print VALUE "Tag_Sequence\tCount\n";
	for (my $j = 0; $j <= $#$degoutput; $j++)
	{
		my $g = $degoutput->[$j];
		if (looks_like_number($g->{"value$i"}))
		{
			print VALUE $g->{GeneNames};
			print VALUE "\t";
			print VALUE $g->{"value$i"}; 
			print VALUE "\n";
		}
	}
	close VALUE;
}

sub raDegoutParse ($)
{
  my ($f) = @_;
  my $line;
  my @deg;
  open DEG, $f or die "cannot open < $f $!";
  $line = <DEG>; chomp $line;
  $line =~ s/"//g;
  my @header = split /\t/, $line;
  while (<DEG>)
  {
    chomp;
    my @e = split /\t/;
    my $g = {};
    for (my $i = 0; $i <= $#header; $i++)
    {
      $g->{$header[$i]} = $e[$i];
    }
    push @deg, $g;
  }
  close DEG;
  return \@deg;
}

__END__

=head1 NAME

edgeR - join two files of degseq and microarray.

=head1 VERSION

edgeR 1.0

=head1 SYNOPSIS

perl map2graph.pl [-degseq file] [-microarray file] 

=head1 DESCRIPTION

edgeR will help you to merge degseq and microarray result.

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
