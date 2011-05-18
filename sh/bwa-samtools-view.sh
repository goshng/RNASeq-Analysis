# Author: Sang Chul Choi
# Date  : 

function bwa-samtools-view {
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
      # bwa samse [-n maxOcc] <in.db.fasta> <in.sai> <in.fq> > <out.sam>
      in_db_fasta=$DATADIR/$GENOMEFASTA-bwa
      in_sai=$DATADIR/SRR031130.sai
      in_fq=$DATADIR/SRR031130.fastq 
      out_sam=$DATADIR/SRR031130.sam
      out_bam=$DATADIR/SRR031130.bam
      $SAMTOOLS view -b -S $out_sam > $out_bam
      echo "Check $DATADIR/SRR031130.bam"

      break
    fi
  done

}
