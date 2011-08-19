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
my $VERSION = 'bwa-sam.pl 1.0';

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
            'rnaseq=s',
            'micro=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $out;
my $outfile;
my $rnaseq;
my $micro;

if (exists $params{out})
{
  $out = $params{out};
  open ($outfile, ">", $out) or die "cannot open > $out: $!";
}
else
{
  $outfile = *STDOUT;   
}

if (exists $params{rnaseq})
{
  $rnaseq = $params{rnaseq};
}
else
{
  &printError("RNA-Seq's result file is missing");
}

if (exists $params{micro})
{
  $micro = $params{micro};
}
else
{
  &printError("Microarray's result file is missing");
}

###############################################################################
# DATA PROCESSING
###############################################################################

if ($cmd eq "deseq")
{
  my %geneResult;
  my $line;
  #############################################################################
  # Read DESeq result
  open DESEQ, $rnaseq or die "cannot open < $rnaseq";
  $line = <DESEQ>;
  while ($line = <DESEQ>)
  {
    my @e = split /\s+/, $line;
    my $genename = $e[1];
    my $log2fc = $e[6];
    my $pval = $e[8];
    my $rec = {};
    $geneResult{$genename} = $rec;
    $rec->{fcrnaseq} = $log2fc;
    $rec->{pvalrnaseq} = $pval;
  }
  close DESEQ;

  #############################################################################
  # Read microarray result
  open MICRO, $micro or die "cannot open < $micro";
  $line = <MICRO>;
  while ($line = <MICRO>)
  {
    my @e = split /\s+/, $line;
    my $genename = $e[4];
    my $log2fc = log($e[3])/log(2);
    my $pval = $e[0];
    unless (exists $geneResult{$genename})
    {
      die "$genename does not exist in microarray";
    }
    my $rec = $geneResult{$genename}; 
    $rec->{fcmicro} = $log2fc;
    $rec->{pvalmicro} = $pval;
  }
  close DESEQ;

  print $outfile "gene\tfcrnaseq\tfcmicro\tpvalrnaseq\tpvalmicro\n";
  foreach my $g (keys %geneResult)
  {
    my $rec = $geneResult{$g};
    if (exists $rec->{pvalmicro})
    {
      print $outfile "$g\t$rec->{fcrnaseq}\t$rec->{fcmicro}\t$rec->{pvalrnaseq}\t$rec->{pvalmicro}\n";
    }
  }
}
if (exists $params{out})
{
  close $outfile;
}
__END__
=head1 NAME

bwa-sam - summary of BWA alignment

=head1 VERSION

bwa-sam 1.0

=head1 SYNOPSIS

perl pl/combine-deseq-microarray.pl deseq -rnaseq output/smutans12/1/bwa/tw1glugal.deseq -micro emails/from/burne/072711/tw1glugal.txt

=head1 DESCRIPTION

Combine DESeq result and microarray results.

  command:
    deseq - rnaseq's result is from DESeq.

=head1 OPTIONS

=over 8

=item B<-rnaseq> <file>

A result file of DESeq.

=item B<-micro> <file>

A file from Burne Lab.

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
