#!/usr/bin/perl

# extractcds
# C. Letondal, letondal@pasteur.fr

use Bio::SeqIO;
use Getopt::Std;
use Getopt::Long;

&init;
&process;

sub usage {
    my $program = `basename $0`;
    chop($program);
    print STDERR "
      $program [ -f format ] [ -p ] [ -l ] [ -g ] [ -h ] seqfile
      Extracts translation(s) from an Embl or Genbank entry

      -f format  : Genbank or Embl format
      -p         : show product
      -g         : show gene
      -l         : show location
      -h         : this message
      seqfile    : entry file (or stdin)
";

}

sub init {
    getopts('f:plgh');
    if ($opt_h) {&usage; exit;}
    if ($opt_f) {$format = $opt_f; } else { $format = "Genbank";}
    if ($opt_p) {$show_product = 1;}
    if ($opt_l) {$show_location = 1;}
    if ($opt_g) {$show_gene = 1;}
    $seqfile = ($ARGV[0]) ? $ARGV[0] : "-";
}

sub process {
    $seqio  = Bio::SeqIO->new (-format => $format , -file => $seqfile);
    $out = Bio::SeqIO->new( '-format' => 'fasta', '-fh' => \*STDOUT );
    while (my $seqobj = $seqio->next_seq())
    {
    foreach $feat ( $seqobj->all_SeqFeatures() ) {
        if ($feat->has_tag('translation')) {
            # if ($feat->has_tag('protein_id')) {
            if ($feat->has_tag('locus_tag')) {
                # @protein_id = $feat->each_tag_value('protein_id');
                @protein_id = $feat->each_tag_value('locus_tag');
                $protein_id = $protein_id[0];
            } elsif ((! $show_gene) && $feat->has_tag('gene')) {
                @gene = $feat->each_tag_value('gene');
                $protein_id = $gene[0];
            } else {
                $protein_id = "";
            }

            if ($show_gene && $feat->has_tag('gene')) {
                @gene = $feat->each_tag_value('gene');
                $gene = $gene[0];
            } else {
                $gene = "";
            }

            if ($show_product && $feat->has_tag('product')) {
                @product = $feat->each_tag_value('product');
                $product = $product[0];
            } else {
                $product = "";
            }

            # $desc = "$protein_id $gene $product";
            $desc = "$protein_id";
            if ($show_location) {
                $desc .= " (" . $feat->start . "," . $feat->end . ")";
            }
            $desc =~ s/\s+/ /g;

            @translation = $feat->each_tag_value('translation');

            my $seq = Bio::Seq->new (
                    -seq  => $translation[0],
                    -id   => $desc
                    );
            $out->write_seq($seq);
        }
    }
    }
}
