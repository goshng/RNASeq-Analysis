# Author: Sang Chul Choi
# Date  : 

function bwa-danko-detecttranscriptsem {
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
      o=$DATADIR/$FUNCNAME.bed

cat>$out_bed_r<<EOF
require(GROseq)
load("$out_bed_rdata")
x <- DetectTranscriptsEM (p=data[,c(1:3,6)])
write.table(x\$Transcripts, "$o", row.names=F, col.names=F, sep="\t", quote=F)
EOF
       
      Rscript $out_bed_r
      echo "Check $o"

      break
    fi
  done

}
