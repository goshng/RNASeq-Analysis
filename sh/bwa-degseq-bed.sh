function bwa-degseq-bed {
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
      READLENGTH=$(grep READLENGTH $SPECIESFILE | cut -d":" -f2)
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {11..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        for i in {1..6}; do
          cut -f$i $DATADIR/$FASTQNUM.bed > $i
        done
        awk '{print "chr1"}' 1 > 11
        awk '{s=$1+"'"$READLENGTH"'";print s}' 2 > 3
        paste 11 2 3 4 5 6 > $DATADIR/$FASTQNUM.bed.degseq
        rm 1 11 2 3 4 5 6
        echo "Check file $DATADIR/$FASTQNUM.bed.degseq"
      done
      break
    fi
  done
}
