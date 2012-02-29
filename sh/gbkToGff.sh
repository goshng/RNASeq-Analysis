#!/bin/bash
#genome=SMU109-OMZ175
#genome=SMU86-U2A
#genome=SMU21-1SM1
#genome=SMU109
#genome=SMU86
genome=SMU21
PERLBIN=/opt/nextgen/bin
rm -rf temp
rm -f $genome.temp.gff
mkdir temp
perl splitgbk.pl raw_data/${genome}_allcontigs.gbf

numcontig=`ls temp | wc -l`

for i in `seq 1 $numcontig`; do
  # genomeSize=`faSize -detailed raw_data/SMU21_6252010.fasta | awk -v contig=$i '$1=="gnl|cornell|contig"contig {print $2}'`
  genomeSize=`cat raw_data/$genome.size | awk -v contig=$i '$1=="gnl|cornell|contig"contig {print $2}'`
  
  perl $PERLBIN/bp_genbank2gff.pl -stdout --file temp/contig$i.gbk |  awk '$1 ~ /^#/ || (NF == 9 && ($4 != 1 || $5 != '$genomeSize')) {print $0}' | sed 's/^unknown/contig'$i'/g' | perl -e '
while (<>) {
    my $line = $_;
    if ($line =~ /^\#/) {
	print $line;
	next;
    }
    chomp($line);
    my @fields=split(); 
    die if (scalar(@fields) != 9);
    my $str=$fields[8]; 
    chomp($str); 
    my $genename;
    my @strs = split(/;/, $str);
    for (my $i=0; $i < scalar(@strs); $i++) {
      if ($strs[$i] =~ /^Parent=([A-Za-z0-9_]+)/) {
        $genename=$1;
      }
      if ($strs[$i] =~ /^ID=([A-Za-z0-9_]+)/) {
        $genename=$1;
      }
      if ($strs[$i] =~ /^locus_tag=([A-Za-z0-9_]+)/) {
        $genename=$1;
      }
    }
   print STDERR "could not find name for $line\n" if (!$genename);
    $fields[8] = $genename;
    print join("\t", @fields) . "\n";
}' >> $genome.temp.gff
done

awk '$3=="gene" {start=$4; end=$5}; $3=="CDS" {if ($4 != start || $5 != end) {print "error "NR}}' $genome.temp.gff
awk 'NR ==1 || ($0 !~ /^#/ && $3=="gene") {print $0}' $genome.temp.gff | sed "s/contig/"$genome".contig/g" >  $genome.gff
awk '$0 !~ /^#/ && $3=="gene" {start=$4; end=$5}; $3=="CDS" {if ($4 != start || $5 != end) {print "error "NR}}' $genome.gff
rm -f $genome.temp.gff
rm -rf temp
