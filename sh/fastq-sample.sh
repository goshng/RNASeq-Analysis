
function fastq-sample {
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
      echo -n "In what intervals do you want to sample? (e.g., 100 for 1%) "
      read INTERVAL

      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        perl pl/$FUNCNAME.pl \
          --fastq $GZIPFASTAQFILE \
          --out $DATADIR/$FASTQNUM.subsample \
          --interval $INTERVAL
        gzip $DATADIR/$FASTQNUM.subsample-$INTERVAL
        echo "Check $DATADIR/$FASTQNUM.subsample-$INTERVAL.gz"
        break
      done

      break
    fi
  done

}
