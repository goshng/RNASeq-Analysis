# Author: Sang Chul Choi
# Date  : 

function bwa-samse {
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
      in.db.fasta=$DATADIR/$GENOMEFASTA-bwa
      in.sai=$DATADIR/SRR031130.sai
      in.fq=$DATADIR/SRR031130.fastq 
      out.sam=$DATADIR/SRR031130.sam
      $BWA samse -n 1 \
        -f $DATADIR/SRR031130.sam \
        $DATADIR/$GENOMEFASTA-bwa \
        $DATADIR/SRR031130.sai \
        $DATADIR/SRR031130.fastq 
      echo "Check $DATADIR/SRR031130.sam"

      break
    fi
  done

}
