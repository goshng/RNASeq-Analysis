# Author: Sang Chul Choi
# Date  : 

function fastxtoolkit-fastq_to_fasta {
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

      in_fq=$DATADIR/SRR03113S.fastq 
      out_fa=$DATADIR/SRR03113S.fasta
      fastq_to_fasta -i $in_fq -o $out_fa
      echo "Check $out_fa"

      break
    fi
  done

}
