#!/opt/local/bin/perl -w
use strict;
use Bio::SeqIO;
  
my $file         = shift; # get the file name, somehow
my $seqio1_object = Bio::SeqIO->new(-file => $file, -format => 'fasta');
$file         = shift; # get the 2nd file name, somehow
my $seqio2_object = Bio::SeqIO->new(-file => $file, -format => 'fasta');

while (my $seq1_object = $seqio1_object->next_seq) {
  my $seq2_object = $seqio2_object->next_seq;
  unless ($seq1_object->seq() eq $seq2_object->seq())
  {
    print $seq1_object->length,"\t",$seq1_object->length,"\n";
  }
}
print "done!\n";
