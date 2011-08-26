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
my $VERSION = 'srnapredict2bed.pl 1.0';

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
            'srnapredict=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $srnapredict;
my $out;
my $outfile;

if (exists $params{srnapredict}) 
{
  $srnapredict = $params{srnapredict};
}
else
{
  &printError("srnapredict file is missing");
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

###############################################################################
# DATA PROCESSING
###############################################################################
if ($cmd eq "bed")
{
  # From 7-th line
  # ************************************************************************************************
  my $line;
  open SRNA, $srnapredict or die "cannot open $srnapredict $!";
  for (my $i = 0; $i < 6; $i++)
  {
    $line = <SRNA>;
  }
  while ($line = <SRNA>)
  {
    chomp $line;
    next if $line =~ /^~/;
    last if length($line) == 0;
    last if $line =~ /^\*+/;
    my @e = split /\t/, $line;
    my $chrom = "chr1";
    my $chromStart = $e[9];
    my $chromEnd = $e[10];
    my $name = $e[0];
    my $score = 0;
    my $strand = "+";
    print $outfile "$chrom\t$chromStart\t$chromEnd\t$name\t$score\t$strand\n"; 
    #print $outfile "$chrom\t$chromStart\t$chromEnd\t$name\t$score\t$strand\t$chromStart\t$chromEnd\t255,0,0\n"; 
  }
  close SRNA;
}

if (exists $params{out})
{
  close $outfile;
}
__END__
=head1 NAME

srnapredict2bed - Subsample a FASTQ file

=head1 VERSION

srnapredict2bed 1.0

=head1 SYNOPSIS

perl srnapredict2bed.pl bed -srnapredict file -out file

=head1 DESCRIPTION

srnapredict2bed will help you to convert sRNAPredict file to a BED file.

=head1 OPTIONS

command: sample, cut
  bed - convert sRNAPredict file to a BED file

=over 8

=item B<-srnapredict> <file>

An sRNAPredict file.

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


