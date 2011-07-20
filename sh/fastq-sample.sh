
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

      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        GZIPFASTAQFILE=$(grep $FASTQNUM $SPECIESFILE | cut -d":" -f2)
        perl pl/$FUNCNAME.pl \
          --fastq $GZIPFASTAQFILE \
          --out $DATADIR/$FASTQNUM.subsample \
          --interval 1000
        gzip $DATADIR/$FASTQNUM.subsample
        echo "Check $DATADIR/$FASTQNUM.subsample.gz"
      done

      break
    fi
  done

}
