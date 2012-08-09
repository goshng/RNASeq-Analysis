# Create tmp/x.tab file.
R --vanilla <<RSCRIPT
library(rtracklayer)
library(GenomicFeatures)
gffFile <- "data/analysis/NC_004350.gff"
sm.gff <- import.gff3(gffFile)
sm.gene <- sm.gff[sm.gff\$type=="gene",]
x <- data.frame(a=sm.gene\$locus_tag,b=sm.gene\$old_locus_tag)
write.table(x,file="/tmp/x.tab",quote=F,col.names=F)
RSCRIPT

# Check if 
perl -e '
while (<>) {
  my @a = split /\s+/;
	my $new = $a[1];
	my $old = $a[2];
	if ($old eq 'NA')
	{
	  print "$new\n"; 
	}
	else
	{
	  $new =~ s/\_/\./;
		if ($new ne $old)
		{
		  print "$old is not found in the new set\n";
		}
	}
}
' < /tmp/x.tab > /tmp/y.tab

# SMU_2166, and SMU_2167 are two new genes and the others are tRNAs or rRNAs.
echo See /tmp/y.tab what new genes are NA in the old set of genes.
