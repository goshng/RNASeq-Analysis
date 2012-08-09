#!/opt/local/bin/perl -w
###############################################################################
# Copyright (C) 2012 Sang Chul Choi
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
use Bio::SeqIO;

$| = 1; # Do not buffer output
my $VERSION = 'bioperl-get-protein.pl 1.0';

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
            'in=s',
            'out=s',
            '<>' => \&process
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

my $out;
my $outfile;
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
my $gbkfile=$params{in};

my $inseq = Bio::SeqIO->new(
                             -file   => "$gbkfile",
                             -format => "Genbank"
                           );

$out = Bio::SeqIO->new( '-format' => 'fasta', '-fh' => $outfile );

while (my $seq = $inseq->next_seq) {
  # print $seq->write_GFF(),"\n";

  foreach my $feat ( $seq->get_SeqFeatures() ) {
    my $locus_tag = "";
    my $translation = "";
    my @a;
    foreach my $tag ( $feat->get_all_tags() ) {
      if ($tag eq 'locus_tag') {
        @a = $feat->get_tag_values($tag);
        $locus_tag = $a[0];
      }
      if ($tag eq 'translation') {
        @a = $feat->get_tag_values($tag);
        $translation = $a[0];
      }
    }

      #print $outfile ">$locus_tag\n$translation\n";
    #}

    if ($locus_tag ne "" and $translation ne "") {
      my $seq = Bio::Seq->new (
                               -seq  => $translation,
                               -id   => $locus_tag
                              );
      $out->write_seq($seq);
    }
  }
}

if (exists $params{out}) {
  close $outfile;
}
__END__
=head1 NAME

bioperl-get-protein - Extract protein sequences from a genbank file

=head1 VERSION

bioperl-get-protein 1.0

=head1 SYNOPSIS

perl pl/bioperl-get-protein.pl -in file.gbk -out file.fa

=head1 DESCRIPTION

bioperl-get-protein will help you to extract protein sequences from a genbank
file.

=head1 OPTIONS

command: sample, cut
  sample - sample short reads
  cut    - truncate tailing nucleotides from short reads
  keep   - keep only short reads with minimum length
  remove - remove reads greater than or equal to a given length
  max    - keep short reads with maximum length

=over 8

=item B<-in> <file>

A genbank file name.

=item B<-out> <number>

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
