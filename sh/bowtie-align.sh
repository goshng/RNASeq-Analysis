function bowtie-align {
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
        COMMAND="zcat $GZIPFASTAQFILE | \
                 BOWTIE_INDEXES=$BOWTIEDIR \
                 $BOWTIE \
                 -t \
                 --solexa1.3-quals \
                 --threads $NUMBERCPU \
                 $GENOMEFASTA \
                 - $BOWTIEDIR/$FASTQNUM.bowtie"

        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND >> $BATCHFILE
        else
          echo $COMMAND >> now
        fi
      done

      break
    fi
  done
}
