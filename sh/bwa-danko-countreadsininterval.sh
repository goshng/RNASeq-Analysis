# Author: Sang Chul Choi
# Date  : 

function bwa-danko-countreadsininterval {
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
      out_bed_danko=$DATADIR/SRR031130.danko
      cut -f1,2,3,6 $out_bed > $out_bed_danko
      echo "Check $out_bed_danko"

      break
    fi
  done

}
