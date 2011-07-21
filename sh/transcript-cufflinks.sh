function transcript-cufflinks {
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

      echo -n "Do you wish to run a batch? (e.g., y/n)"
      read WISH
      if [ "$WISH" == "y" ]; then
        BATCH=YES
        BATCHFILE=batch.sh
        echo "#!/bin/bash" > $BATCHFILE
      fi

      GENOMEFASTA=$(basename $REFGENOMEFASTA)
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        COMMAND1="BOWTIE_INDEXES=$BOWTIEDIR/ \
                  tophat --num-threads 4 --solexa1.3-quals \
                  -o $BOWTIEDIR/tophat$FASTQNUM \
                  $GENOMEFASTA $GZIPFASTAQFILE"
        COMMAND2="cufflinks -o $BOWTIEDIR/cufflinks$FASTQNUM \
                  $BOWTIEDIR/tophat$FASTQNUM/accepted_hits.bam"
        # grep "RefSeq use ctrl+v and tab to insert a tab gene" NC_004350.gff > 1.gff
        # Many duplicates are manually edited.
        # 1.gff must be made. See data/1.gff file.
        COMMAND3="cuffcompare -r 1.gff \
                  -o $BOWTIEDIR/cuffcompare$FASTQNUM \
                  $BOWTIEDIR/cufflinks$FASTQNUM/transcripts.gtf"

        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND1 >> $BATCHFILE
          echo $COMMAND2 >> $BATCHFILE
          # echo $COMMAND3 >> $BATCHFILE
        else
          echo $COMMAND1 | bash
          echo $COMMAND2 | bash
          # echo $COMMAND3 | bash
          echo "Check $BOWTIEDIR/tophat"
          echo "Check $BOWTIEDIR/cufflinks"
        fi
        break
      done
      break
    fi
  done
}
