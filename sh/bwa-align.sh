# Author: Sang Chul Choi
# Date  : 

function bwa-align {
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
      $BWA aln -t $NUMBERCPU $DATADIR/$GENOMEFASTA-bwa $DATADIR/SRR031130.fastq > $DATADIR/SRR031130.sai
      echo "Check $DATADIR/SRR031130.sai"

      break
    fi
  done

}
