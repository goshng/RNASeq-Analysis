#!/usr/bin/perl
#	Program outputs the column specificed from a TSV file (0 = first column ...)
#
#
#	Assume that repeat reads are already stripped -- can do this by piping to grep -v 
#	e.g.: samtools view INPUT.bam | grep "XA:Z" | sam2bed
#	That removes the 2nd best match.
#
#	Alternatively, threshold the mapping quality to be >0 (this is done below).
#
# 1  QNAME   String 
# 2  FLAG    Int 
# 3  RNAME   String 
# 4  POS     Int
# 5  MAPQ    Int 
# 6  CIGAR   String 
# 7  RNEXT   String 
# 8  PNEXT   Int 
# 9  TLEN    Int
# 10 SEQ     String 
# 11 QUAL    String
#
# An input file looks like this:
# HWI-ST397_0000:3:2208:11124:55733#CGATGT        0
# gi|24378532|ref|NC_004350.1|    20      37      100M    *       0       0
# CACTTTTCCACAAGAAAAGATGCTATCGAATCTCTTGATTAACAGAATTTATCATCTTTTCCACAAATTGTGGAAAACTTATGTCCACATTGTGGACTCA
# HGGHGHHHHHHHGHHHHGGHGHGHHHGHHHHHHHFHHFHFHHHFFDEHHH2DDDDDHFHHHHABHHEEED9C@CBGFFGEHHFDHEFFEHFBHHEEDBC?
# XT:A:U  NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:100

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
