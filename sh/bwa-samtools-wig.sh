function bwa-samtools-wig {
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

      REFGENOMELENGTH=$(grep REFGENOMELENGTH $SPECIESFILE | cut -d":" -f2)
      READLENGTH=$(grep READLENGTH $SPECIESFILE | cut -d":" -f2)
      NUMFASTQFILE=$(grep NUMFASTQFILE $SPECIESFILE | cut -d":" -f2)
      for g in $(eval echo {1..$NUMFASTQFILE}); do
        FASTQNUM=FASTQ$(printf "%02d" $g)
        perl pl/$FUNCNAME.pl \
          -genomeLength $REFGENOMELENGTH \
          -readLength $READLENGTH \
          -bed $DATADIR/$FASTQNUM.bed \
          -wig $DATADIR/$FASTQNUM.wig
      done
      break
    fi
  done
}

