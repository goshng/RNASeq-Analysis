# Author: Sang Chul Choi
# Date  : 

function bwa-danko-metagene {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION
      global-variable $SPECIES $REPETITION
      read-species

      GENOMEFASTA=$(basename $REFGENOMEFASTA)
      in_db_fasta=$DATADIR/$GENOMEFASTA-bwa
      in_sai=$DATADIR/SRR031130.sai
      in_fq=$DATADIR/SRR031130.fastq 
      out_sam=$DATADIR/SRR031130.sam
      out_bam=$DATADIR/SRR031130.bam
      out_bam_sorted=$DATADIR/SRR031130.sorted
      out_bed=$DATADIR/SRR031130.bed
      out_bed_r=$DATADIR/SRR031130-$FUNCNAME.R
      out_bed_rdata=$DATADIR/SRR031130.RData
      out_bed_danko=$DATADIR/SRR031130.danko
      f=$DATADIR/feature-genome.out-genestart
      o=$DATADIR/$FUNCNAME.out
      ops=$DATADIR/$FUNCNAME.out.ps

cat>$out_bed_r<<EOF
require(GROseq)
load("$out_bed_rdata")
y <- read.table("$f")
x <- MetaGene (p=data[,c(1:3,6)], f=y, size=3, up=5000)
write.table(x, "$o")
# plot(as.integer(row.names(x))-5000,x\$x, xlim=c(-1000,1000))
x <- read.table("$o")
postscript ("$ops",  width=10, height=10, horizontal = FALSE, onefile = FALSE, paper = "special")
plot(as.integer(row.names(x))-5000,x\$x, xlim=c(-1000,1000),cex=0.5,pch=21)
dev.off()
EOF
       
      Rscript $out_bed_r
      echo "Check $out_bed_r and $o"
      echo "open $ops"

      break
    fi
  done

}
