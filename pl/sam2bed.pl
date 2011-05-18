#!/usr/bin/perl
#	Program outputs the column specificed from a TSV file (0 = first column ...)
#
#
#	Assume that repeat reads are already stripped -- can do this by piping to grep -v 
#	e.g.: samtools view INPUT.bam | grep "XA:Z" | sam2bed
#	That removes the 2nd best match.
#
#	Alternatively, threshold the mapping quality to be >0 (this is done below).

while(<STDIN>) {
	@SPL = split(/[\s\t]/);
	$str  = $SPL[1];
	$chrom = $SPL[2];
	$chromStart = $SPL[3];
	$mapq = $SPL[4];
	if($str == 0 && $mapq>0) {
	  print "$chrom\t$chromStart\t$chromStart\tn\t$mapq\t+\n"; 
	}
	if($str == 16 && $mapq>0) {
          print "$chrom\t$chromStart\t$chromStart\tn\t$mapq\t-\n";
	}
}
