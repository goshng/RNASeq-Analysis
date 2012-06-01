###############################################################################
# Copyright (C) 2012 Sang Chul Choi
#
# This file is part of RNASeq Analysis.
# 
# RNASeq Analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RNASeq Analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RNASeq Analysis.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

function convert-gff2txdb {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable

      REFGENOMEGFF=$(grep ^REFGENOMEGFF\: $SPECIESFILE | cut -d":" -f2)
      REFGENOMETXDB=$(grep ^REFGENOMETXDB\: $SPECIESFILE | cut -d":" -f2)
      REFGENOMETXDBBASE=$(basename $REFGENOMETXDB)

      RTEMP=$BWADIR/$RANDOM.R
      COMMAND="Rscript $RTEMP"
cat>$RTEMP<<EOF
library(DESeq)
library(ShortRead)
library(rtracklayer)
library(GenomicRanges)
library(VariantAnnotation)
library(GenomicFeatures)
gffFile <- "$REFGENOMEGFF"
sm.gff <- import.gff3(gffFile)
sm.source <- sm.gff[1,] #sm.gff$type=="region",]
# For AEXT01 contig genome
# sm.source <- sm.gff[grep("^AEXT01",sm.gff$ID),] 
sm.gene <- sm.gff[sm.gff\$type=="gene",]
sm.exon <- sm.gff[sm.gff\$type=="exon",]
sm.CDS <- sm.gff[sm.gff\$type=="CDS",]
# For AEXT01 contig genome
sm.CDS <- sm.gff[sm.gff\$type=="mRNA",]
sm.rrna <- sm.gff[sm.gff\$type=="rRNA",]

transcripts <- 
  data.frame( tx_id=seq(length(sm.gene\$locus_tag)),
              tx_name=sm.gene\$locus_tag,
              tx_chrom=sm.gene\$space, 
              tx_strand=sm.gene\$strand,
              tx_start=start(sm.gene),
              tx_end=end(sm.gene) )

# Compare CDS and gene locus tags
x <- c()
y <- c()
z <- unlist(sm.CDS\$Parent)
for (i in sm.gene\$ID) {
  if (sum(z == i) > 0) {
    if (length(start(sm.CDS)[z == i]) != 1) {               
      print(paste("Check",gffFile))
      print(start(sm.CDS)[z == i])                          
      stop(paste("There are multiple CDS for",i))                               
    }
    x <- c(x,start(sm.CDS)[z == i])
    y <- c(y,end(sm.CDS)[z == i])
  } else {
    x <- c(x,NA)
    y <- c(y,NA)
  }
}
grCDS.start <- x
grCDS.end <- y
rm(i,x,y,z)

# 
splicings <- 
  data.frame( tx_id=seq(length(sm.gene\$locus_tag)),
              exon_rank=seq(length(sm.gene\$locus_tag)),
              exon_start=start(sm.gene),
              exon_end=end(sm.gene),
              exon_name=sm.gene\$locus_tag,
              cds_start=grCDS.start,
              cds_end=grCDS.end )

chrominfo <-
  data.frame( chrom=names(sm.gene),
              length=end(sm.source),
              is_circular=FALSE )

TxDb.Smutans.UA159.uid57947.knownGene <- makeTranscriptDb (transcripts, splicings, chrominfo=chrominfo)
txdb <- TxDb.Smutans.UA159.uid57947.knownGene 
saveFeatures(txdb,file="$REFGENOMETXDB")
print("Check $REFGENOMETXDB")

EOF
      if [ "$BATCH" == "YES" ]; then
        echo $COMMAND >> $BATCHFILE
      else
        #echo $COMMAND | bash
        #rm $RTEMP
        #echo $RTEMP
        echo "Edit and run $RTEMP"
	echo "e.g., Rscript $RTEMP"
      fi
      break
    fi
  done
}

function convert-bed2txdb {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable

      REFGENOMEBED=$(grep ^REFGENOMEBED\: $SPECIESFILE | cut -d":" -f2)
      if [ "$REFGENOMEBED" == "" ]; then 
        echo No REFGENOMEBED in $SPECIESFILE; break; 
      fi
      REFGENOMETXDB=$(grep ^REFGENOMETXDB\: $SPECIESFILE | cut -d":" -f2)
      if [ "$REFGENOMETXDB" == "" ]; then 
        echo No REFGENOMETXDB in $SPECIESFILE; break; 
      fi

      RTEMP=$BWADIR/$RANDOM.R
      COMMAND="Rscript $RTEMP"
cat>$RTEMP<<EOF
library(GenomicFeatures)
bedFile <- "$REFGENOMEBED"
sm.bed <- read.table(bedFile,header=F)
sm.bed\$V2 <- sm.bed\$V2 + 1
transcripts <- 
  data.frame( tx_id=seq(length(sm.bed\$V1)),
              tx_name=sm.bed\$V4,
              tx_chrom=sm.bed\$V1,
              tx_strand='+',
              tx_start=sm.bed\$V2,
              tx_end=sm.bed\$V3 )

splicings <- 
  data.frame( tx_id=seq(length(sm.bed\$V1)),
              exon_rank=seq(length(sm.bed\$V1)),
              exon_start=sm.bed\$V2,
              exon_end=sm.bed\$V3,
              exon_name=sm.bed\$V4,
              cds_start=sm.bed\$V2,
              cds_end=sm.bed\$V3 )

chrominfo <-
  data.frame( chrom="chr1",
              length=sm.bed[length(sm.bed\$V1),]\$V3,
              is_circular=FALSE )

txdb <- makeTranscriptDb (transcripts, splicings, chrominfo=chrominfo)
saveFeatures(txdb,file="$REFGENOMETXDB")
print("Check $REFGENOMETXDB")
EOF
      if [ "$BATCH" == "YES" ]; then
        echo $COMMAND >> $BATCHFILE
      else
        echo $COMMAND | bash
        rm $RTEMP
        echo $RTEMP
        echo "Edit and run $RTEMP"
	echo "e.g., Rscript $RTEMP"
      fi
      break
    fi
  done
}
