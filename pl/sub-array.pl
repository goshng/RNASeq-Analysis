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

sub arrayNumberSlope (@) {
  my @a = @_;
  my $v = 0;
  my $v1 = $a[0];
  my $d = 0;
  for (my $i = 1; $i <= $#a; $i++)
  {
    $d++;
    my $v2 = $a[$i];
    if ($v2 > 0)
    {
      $v += (($v1 - $v2)/$d);
      $v1 = $v2;
      $d = 0;
    }
  }
  if ($v < 0)
  {
    $v = 0;
  }
  return $v;
}

sub arrayNumberNegate (@) {
  my @a = @_;
  my $v = 0;
  for (my $i = 0; $i <= $#a; $i++)
  {
    $a[$i] = -$a[$i]
  }
  return @a;
}

sub arrayNumberNonzero (@) {
  my @a = @_;
  my $v = 0;
  for (my $i = 0; $i <= $#a; $i++)
  {
    if ($a[$i] != 0)
    {
      $v++;
    }
  }
  return $v;
}

1;
