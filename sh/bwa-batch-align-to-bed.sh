# Author: Sang Chul Choi

function bwa-batch-align-to-bed {
  PS3="Choose the species for $FUNCNAME: "
  select SPECIES in ${SPECIESS[@]}; do 
    if [ "$SPECIES" == "" ];  then
      echo -e "You need to enter something\n"
      continue
    else  
      echo -n "What repetition do you wish to run? (e.g., 1) "
      read REPETITION

      BATCH=YES
      BATCHFILE=batch.sh
      echo "#!/bin/bash" > $BATCHFILE
      bwa-align
      bwa-samse
      bwa-samtools-view
      bwa-samtools-sort
      bwa-samtools-bed

      break
    fi
  done

}
