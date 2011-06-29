# Author: Sang Chul Choi

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
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        COMMAND="$BWA aln -I -t $NUMBERCPU \
                 $DATADIR/$GENOMEFASTA-bwa \
                 $GZIPFASTAQFILE > $DATADIR/$FASTQNUM.sai"

        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND >> $BATCHFILE
        else
          $COMMAND
          echo "Check $DATADIR/$FASTQNUM.sai"
        fi

      done

      break
    fi
  done

}
