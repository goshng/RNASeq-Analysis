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

sub rnaseqPttParse ($) {
  my ($f) = @_;
  my @v;
  # Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
  open PTT, $f or die "cannot open $f $!";
  my $line = <PTT>;
  $line = <PTT>;
  $line = <PTT>;
  while ($line = <PTT>)
  {
    # print $outfile "chr1\t$t->{start}\t$t->{end}\t$t->{name}\t$t->{score}\t$t->{strand}\n";
    chomp $line;
    my @e = split /\t/, $line;
    my $rec = {};
    # Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
    $rec->{Location} = $e[0];
    $rec->{Strand}   = $e[1];
    $rec->{Length}   = $e[2];
    $rec->{PID}      = $e[3];
    $rec->{Gene}     = $e[4];
    $rec->{Synonym}  = $e[5];
    $rec->{Code}     = $e[6];
    $rec->{COG}      = $e[7];
    $rec->{Product}  = $e[8];
    push @v, $rec;
  }
  close PTT;
  return @v;
}

# Use PTT file to extract the genome length.
sub rnaseqPttGetGenomeLength ($) {
  my ($f) = @_;
  my $v;
  # Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
  open PTT, $f or die "cannot open $f $!";
  my $line = <PTT>;
  $line =~ /(\d+)\.\.(\d+)/;
  $v = $2;
  return $v;
}

1;
