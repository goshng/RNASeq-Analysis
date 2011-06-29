# Author: Sang Chul Choi
# Date  : 

function bwa-samtools-bed {
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

      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        COMMAND="$SAMTOOLS view $DATADIR/$FASTQNUM.sorted.bam \
          | perl pl/sam2bed.pl > $DATADIR/$FASTQNUM.bed"

        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND >> $BATCHFILE
        else
          $COMMAND
          echo "Check $DATADIR/$FASTQNUM.bed"
        fi
      done

      break
    fi
  done

}
