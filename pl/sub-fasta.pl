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
=head1 NAME

sub-fasta.pl - Parsers of FASTA files

=head1 SYNOPSIS

  my $fnaSequence = maFastaParse ($fna);

=head1 VERSION

v1.0, Sat Jul 16 07:20:14 EDT 2011

=head1 DESCRIPTION

FASTA file parser.

=head1 FUNCTIONS

=over 4

=item sub maFastaParse ($)

  Argument 1: FASTA file
  Return: String

=item sub peachFastaLength ($)

  Argument 1: FASTA file
  Return: Sequence name and lengths in a hash

=back

=head1 AUTHOR

Sang Chul Choi, C<< <goshng_at_yahoo_dot_co_dot_kr> >>

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

sub maFastaParse ($)
{
  my ($f) = @_;
  my $line;
  my $sequence;
  open FASTA, $f or die "cannot open < $f: $!";
  $line = <FASTA>;
  unless ($line =~ /^>/)
  {
    die "The first character of the first line must be '>': $line";
  }
  while ($line = <FASTA>)
  {
    chomp $line;
    if ($line =~ /^>/)
    {
      die "Only one sequence in the file";
    }
    $sequence .= $line;
  }
  close FASTA;
  return $sequence;
}

sub peachFastaLength ($)
{
  my ($f) = @_;
  my $line;
  my $sequence = "";
  my %fl;
 
  my $curSequenceName = "";
  open FASTA, $f or die "cannot open < $f: $!";
  while ($line = <FASTA>)
  {
    chomp $line;
    if ($line =~ /^>(\S+)/)
    {
      if (length ($curSequenceName) > 0)
      {
        if (exists $fl{$curSequenceName} and $fl{$curSequenceName} == 0)
	{
          $fl{$curSequenceName} = length ($sequence);
	}
	else
	{
	  die "curSequenceName does not exists: see the $f for validity of FASTA file";
	}
      }
      $curSequenceName = $1;
      $fl{$1} = 0;
      $sequence = "";
    }
    else
    {
      $sequence .= $line;
    }
  }
  $fl{$curSequenceName} = length ($sequence);

  close FASTA;
  return %fl;
}

sub rnaseqFastaCreate ($$$)
{
  my ($f, $header, $seq) = @_;
  $seq =~ s/-//g;
  open OUT, ">", $f or die "cannot open $f $!";
  print OUT ">$header\n";
  print OUT "$seq\n";
  close OUT;
}

sub rnaseqFastaAdd ($$$)
{
  my ($f, $header, $seq) = @_;
  $seq =~ s/-//g;
  open OUT, ">>", $f or die "cannot open $f $!";
  print OUT ">$header\n";
  print OUT "$seq\n";
  close OUT;
}

# Return 1 if there is a sequence with many gaps in the 
# FASTA-format file.
sub rnaseqFastaManyGap ($)
{
  my ($f) = @_;
  my $v = 0;
  my $seq = "";
  my $line;
  open OUT, $f or die "cannot < open $f $!";
  while ($line = <OUT>)
  {
    chomp $line;
    if ($line =~ /^>/)
    {
      if (length $seq > 0)
      {
        # Check the sequence has consecutive gaps.
        my $prevChar = ".";
        my $c = 0;
        my $maxc = 0;
        while ($seq =~ /(.)/g)
        {
          if ($1 eq "-" and $prevChar eq "-")
          {
            $c++;
          }
          else
          {
            $c = 0;
          }
          if ($maxc < $c)
          {
            $maxc = $c;
          }
          $prevChar = $1;
        }
        if ($maxc > 40)
        {
          $v = 1;
          last;
        }
      }
      $seq = "";
    }
    else
    {
      $seq .= $line;
    }
  }
  close OUT;
  return $v;
}

1;
