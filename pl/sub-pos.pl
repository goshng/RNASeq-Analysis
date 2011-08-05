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

sub rnaseqPosParse ($) {
  my ($ingene) = @_;
  my @genes;
  my $s = 0;
  open INGENE, "$ingene" or die "cannot open < $ingene";
  while (<INGENE>)
  {
    chomp;
    my @e = split /\t/;
    my $rec = {};
    $rec->{name}  = $e[0]; 
    $rec->{chr}   = $e[1]; 
    $rec->{start} = $e[2]; 
    $rec->{end}   = $e[3]; 
    $rec->{count} = 0;
    unless ($rec->{start} < $rec->{end})
    {
      die "feature start must be less than its end";
    }
    push @genes, $rec;
    $s++;
    if ($s % 100000 == 0)
    {
      print STDERR "$s\r";
    }
  }
  close INGENE;
  return @genes;
}
1;
