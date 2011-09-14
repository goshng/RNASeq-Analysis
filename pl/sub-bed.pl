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

sub rnaseqBedParse ($) {
  my ($f) = @_;
  my @v;
  open BED, $f or die "cannot open $f $!";
  while (<BED>)
  {
    # print $outfile "chr1\t$t->{start}\t$t->{end}\t$t->{name}\t$t->{score}\t$t->{strand}\n";
    chomp;
    my @e = split /\t/;
    my $rec = {};
    $rec->{chr}    = $e[0];
    $rec->{start}  = $e[1];
    $rec->{end}    = $e[2];
    $rec->{name}   = $e[3];
    $rec->{score}  = $e[4];
    $rec->{strand} = $e[5];
    push @v, $rec;
  }
  close BED;
  return @v;
}

1;
