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

      # bwa samse [-n maxOcc] <in.db.fasta> <in.sai> <in.fq> > <out.sam>
      GENOMEFASTA=$(basename $REFGENOMEFASTA)
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        COMMAND="$BWA samse -n 1 \
          -f $BWADIR/$FASTQNUM.sam \
          $BWADIR/$GENOMEFASTA-bwa \
          $BWADIR/$FASTQNUM.sai \
          $GZIPFASTAQFILE"
        if [ "$BATCH" == "YES" ]; then
          echo $COMMAND >> $BATCHFILE
        else
          echo $COMMAND | bash
          echo "Check $BWADIR/$FASTQNUM.sam"
        fi
      done

      break
    fi
  done

}
